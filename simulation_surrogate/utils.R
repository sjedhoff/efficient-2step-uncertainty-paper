################################################################################
##                                                                            ##
##                    Utils for logistic and PCE surrogate                    ##
##                    for the T-Step, I-step, PSIS + IWMM                     ##
##                                                                            ##
################################################################################

library(cmdstanr)
library(here)
library(orthopolynom)
library(loo)

##  T Step PCE:
#   Input: 
#     - config: list
#     - fct: function of the true simulator, dependend on x1
#   Output: list
#     - tmodel_coeffs: matrix with nrow = config$itersampling * config$chains
#     - pce_vars: list with l_poly_coeffs_mat, comb
t_step_pce <- function(config, fct){
  Rcpp::sourceCpp("PCE_helpers.cpp")
  
  ## Get T-Data
  get_pce_vars <- function(p, M, idx=NULL){
    l_polys <- legendre.polynomials(p, normalized=TRUE)
    l_poly_coeffs <- polynomial.coefficients(l_polys)
    l_poly_coeffs_mat <- matrix(0, p+1, p+1)
    for (i in 1:(p+1)){
      for (j in 1:length(l_poly_coeffs[[i]])){
        l_poly_coeffs_mat[i, j] = l_poly_coeffs[[i]][j]
      }
      if (i+1<=p+1){
        l_poly_coeffs_mat[i, seq(i+1, p+1)] <- 0
      }
    }
    # Comb & Poly selection
    comb <- poly_idx_cpp(p, M)
    # first column is the constant polynomial f0
    comb <- comb[-1, , drop = FALSE]
    rownames(comb) <- seq_len(nrow(comb))
    if (!is.null(idx)) {
      # select only desired polynomials
      comb <- comb[idx, , drop = FALSE]
    }
    list("l_poly_coeffs_mat" = l_poly_coeffs_mat,
         "comb" = comb)
  }
  
  pce_vars <- get_pce_vars(p = config$poly_degree, M = config$M, idx = NULL)
  l_poly_coeffs_mat <- pce_vars[[1]]
  comb <- pce_vars[[2]]
  
  w_sim <- seq(-1, 1, length.out=config$n_sims)

  tdata <- list(
    N_sim = config$n_sims,
    sigma_sim_lower = 0.0,
    d = config$poly_degree,
    M = config$M,
    w_sim = matrix(c(w_sim), ncol=config$M),
    y_sim = fct(w_sim) + rnorm(n = config$n_sims, mean = 0, sd = config$sigma_sim),
    l_poly_coeffs = t(l_poly_coeffs_mat),
    comb = comb,
    N_comb = nrow(comb),
    prior_only = 0
  )
  
  ## Fit the T-Model
  tmodel_file <- "stan_code/multi_legendre_pce.stan"
  tmodel <- cmdstan_model(tmodel_file)
  
  tmodel_fit <- tmodel$sample(
    data = tdata,
    seed = config$seed,
    chains = config$chains,
    parallel_chains = config$chains,
    iter_sampling = config$S / config$chains,
    adapt_delta = config$adapt_delta,
    refresh=0
  )
  
  tmodel_coeffs <- tmodel_fit$draws(c("c_0", "c", "sigma_sim"), format = "matrix")
  
  return(list("tmodel_coeffs" = tmodel_coeffs,
              "pce_vars" = pce_vars))
}

#  T Step Logistic:
#   Input: 
#     - config: list
#     - fct: function of the true simulator, dependend on x1
#   Output: list
#     - tmodel_coeffs: matrix with nrow = config$itersampling * config$chains
#     - pce_vars: NULL, just for consistency with the PCE t step
t_step_logistic <- function(config, fct){
  ## The input data from the model M: D_T
  #   - get n_sims datapoints w_sim between -1 and 1
  #   - get the corresponding y values (put in the w in the logistic function) + error
  w_sim <- seq(-1, 1, length.out=config$n_sims)
  tdata <- list(
    N_sim = config$n_sims,
    y_sim = fct(w_sim) + rnorm(n = config$n_sims, mean = 0, sd = config$sigma_sim),
    w_sim = matrix(w_sim, ncol = 1),
    sigma_sim_lower = 0,
    prior_only = 0
  )
  
  ## Fit the T-Model
  stan_file_t <- paste0("stan_code/true_model_logistic_4params.stan")
  tmodel <- cmdstan_model(stan_file_t)
  
  tmodel_fit <- tmodel$sample(
    data = tdata,
    seed = config$seed,
    chains = config$chains,
    parallel_chains = config$chains,
    iter_sampling = config$S / config$chains,
    adapt_delta = config$adapt_delta,
    refresh=0 
  )
  
  tmodel_coeffs <- tmodel_fit$draws(c("c_0", "c", "sigma_sim"), format = "matrix")
  return(list("tmodel_coeffs" = tmodel_coeffs,
              "pce_vars" = NULL))
}


##  Get I-Data
#   Input: 
#     - config
#     - fct: function of the true simulator, dependend on x1
#   Output: 
#     - idata: list with w_exp_1 and y_exp
get_idata <- function(config, fct){
  w_exp <- rep(config$w_exp_gts, config$n_exp)
  idata <- list(
    w_exp_1 = matrix(w_exp, ncol = config$n_exp),
    y_exp = matrix(fct(w_exp)+rnorm(n=config$n_exp, mean=0, sd=config$sigma_exp),
                   ncol = config$n_exp)
  )
  
  return(idata)
}

## Select T-draws for fitting the I-Model:
#  Input: 
#     - config
#     - tmodel_coeffs: T-draws matrix(S, config$poly_degree + 2)
#     - idata
#     - method: character, "PCE" or "logistic"
# Output: 
#     - t_draws_proposal_id: numeric vector with selected ids
#     - t_draws_proposal: tmodel_coeffs matrix with selected rows
#     - imodel_input: list: ids, tmodel_coeffs
select_t_draws <- function(config, tmodel_coeffs, idata, method = "PCE"){
  
  
  w_draws <- truncnorm::rtruncnorm(n = config$n_prior, a = -1, b = 1, mean = 0, sd = 0.5)
  
  ll <- matrix(NA, ncol = config$n_prior, nrow = nrow(tmodel_coeffs))
  
  ## PCE
  if(method == "PCE"){
    l_polys <- legendre.polynomials(config$poly_degree)
    
    for(s in 1:nrow(tmodel_coeffs)){
      for(i in 1:config$n_prior){
        w <- w_draws[i]
        x <- sapply(1:(config$poly_degree+1), function(i) l_polys[[i]] %*% w^(0:(i-1)) )
        mean <- sum(tmodel_coeffs[s,1:(config$poly_degree + 1)] * x)
        sd <- as.numeric(tmodel_coeffs[s,(config$poly_degree + 2)])
        ll[s,i] <- sum(dnorm(idata$y_exp, mean = mean, sd = sd, log = TRUE))
      }
    }
  }
  
  ## Logistic
  else{
    if(method == "logistic"){
      logistic_fct_mean_vec <- function(c_0, c, w){
        return( c_0 / (1 + exp(-c[1] * (w - c[2]))) + c[3])
      }
      
      for(s in 1:nrow(tmodel_coeffs)){
        for(i in 1:config$n_prior){
          mean <- logistic_fct_mean_vec(tmodel_coeffs[s,1], tmodel_coeffs[s,2:4], w_draws[i])
          sd <- tmodel_coeffs[s,5]
          ll[s,i] <- sum(dnorm(idata$y_exp, mean = mean, sd = sd, log = TRUE))
        }
      }
      
    }
    else{
      stop("method must either be 'PCE' or 'logistic'")
    }
  }
  
  ll_means <- rowMeans(ll)
  
  if(config$M_mix == 1){
    probs <- 0.5
  }
  else{
    probs <- seq(0,1,length.out = config$M_mix)
  }
  t_draws_proposal_id <-  which(ll_means %in% quantile(ll_means, probs = probs, type = 1))
  

  t_draws_proposal <- tmodel_coeffs[t_draws_proposal_id,]
  
  imodel_input <- list("ids" = as.numeric(rownames(tmodel_coeffs[t_draws_proposal_id,])),
                        "tmodel_coeffs" = t_draws_proposal)
  
  return(list("t_draws_proposal_id" = as.numeric(rownames(tmodel_coeffs[t_draws_proposal_id,])),
              "t_draws_proposal" = t_draws_proposal,
              "imodel_input" = imodel_input))
}

## Fit I-Model
#  Input:
#     - input: list: ids numeric(N_mixture), tmodel_coeffs matrix(N_mixture, config$poly_degree+2)
#     - imodel: cmdstan_model
#     - pce_vars: list with l_poly_coeffs_mat, comb, if NULL, instead of PCE the logistic one will be used
#     - idata
#     - config
#     - ll_weight: numeric, used to weigh log-likelihoods
# Output: 
#     - w_exp, sigma_exp, log_lik_proposal: matrix(K, N_mixture)
#     - log_lik_evals: number of evaluations of the log probability evaluations
#     - input: (see Input)
#     - N_mixture: numeric, number of T-draws fitted
#     - K: number of T-posterior draws
#     - fits: list of fitted cmdstanr models
fit_imodel <- function(input, imodel, pce_vars = NULL, idata, config, 
                       ll_weight = 1){
  fits <- list()
  
  ## PCE
  if(!is.null(pce_vars)){
    input_list <- lapply(1:length(input$ids), function(i)     imodel_input <- list(
      N_exp = nrow(idata$w_exp_1),
      N_measures = ncol(idata$w_exp_1),
      M = config$M,
      y_exp = idata$y_exp,
      d = config$poly_degree,
      
      l_poly_coeffs = t(pce_vars$l_poly_coeffs_mat),
      N_comb = nrow(pce_vars$comb),
      comb = pce_vars$comb,
      c_0 = as.numeric(input$tmodel_coeffs[i,1]),
      c =  as.numeric(input$tmodel_coeffs[i,2:(nrow(pce_vars$comb)+1)]),
      prior_only = FALSE,
      
      ll_weight = ll_weight
      
    ))
  }
  ## Logistic
  else{
    input_list <- lapply(1:length(input$ids), function(i)     imodel_input <- list(
      N_exp = nrow(idata$w_exp_1),
      N_measures = ncol(idata$w_exp_1),
      y_exp = idata$y_exp,

      c_0 = as.numeric(input$tmodel_coeffs[i,1]),
      c =  as.numeric(input$tmodel_coeffs[i,2:4]),
      prior_only = FALSE
      
    ))
  }
  
  ## Fit the i-model for all cluster centers:
  for(i in 1:length(input$ids)) {
    # fit the i-model
    imodel_fit <- imodel$sample(
      data = input_list[[i]],
      seed = config$seed,
      chains = config$chains,
      parallel_chains = config$chains,
      iter_sampling = config$K / config$chains,
      iter_warmup = config$iter_warmup,
      refresh = config$refresh_cmdstan,
      adapt_delta = config$adapt_delta
    )
    
    fits[[i]] <- imodel_fit
  }
  
  
  w_exp <- matrix(unlist(lapply(fits, function(x)
    as.numeric(x$draws("w_exp[1,1]")))),
    ncol = length(input$ids))
  sigma_exp <- matrix(unlist(lapply(fits, function(x)
    as.numeric(x$draws("sigma_exp")))),
    ncol = length(input$ids))
  log_lik_proposal <- matrix(unlist(lapply(fits,function(x)
    as.numeric(x$draws("log_lik")))),
    ncol = length(input$ids))
  
  # Number of log-lik evaluations
  log_lik_evals <- sum(unlist(lapply(fits, function(fit) sum(fit$sampler_diagnostics()))))
  
  return(list("w_exp" = w_exp, 
              "sigma_exp" = sigma_exp, 
              "log_lik_proposal" = log_lik_proposal, 
              "log_lik_evals" = log_lik_evals,
              "input" = input,
              "N_mixture" = length(input$ids),
              "K" = nrow(w_exp),
              "fits" = fits))
}

## Get Proposal values
# Input:
#   - w_exp, sigma_exp: matrix(K, N_measures)
# Output:
#   - list: w_values, sigma_values vector(K)
get_proposal_values <- function(w_exp, sigma_exp){
  
  N_mixture <- ncol(w_exp)
  K <- nrow(w_exp)

  calculate_counts <- function(K, N_mixture) {
    counts <- rep(0, N_mixture)  
    
    count_per_group <- floor(K / N_mixture)
    remainder <- K %% N_mixture
    
    counts[1:remainder] <- count_per_group + 1
    counts[(remainder + 1):N_mixture] <- count_per_group
    
    return(counts)
  }
  
  num_draws <- calculate_counts(K, N_mixture)
  
  proposal_values <- unlist(lapply(1:N_mixture, function(i) sample(w_exp[,i], num_draws[i])))
  proposal_sigma_values <- unlist(lapply(1:N_mixture, function(i) sample(sigma_exp[,i], num_draws[i])))
  
  
  return(list("w_values" = proposal_values,
              "sigma_values" = proposal_sigma_values))
  
}

## Calculation of log ratios using stan:
#  Input:
#     - config
#     - idata: list with w_exp_1 and y_exp
#     - pce_vars: list with l_poly_coeffs_mat and comb, if NULL, the logistic version is used
#     - imodel_input: list ids, tmodel_coeffs (from function select_t_draws_pce)
#     - tmodel_coeffs: matrix(S, config$poly_degree + 2)
#     - proposals: list with w_values and sigma_values (from function get_proposal_values)
#     - K: number of i-draws from the proposal for each t-draw
#     - N_mixture: number of t-draws used for the mixture proposal
# Output:
#     - log_r: log ratios matrix(K, S) 
#     - proposals: list wirh w_values, sigma_values
#     - log_ratio_evals: number of log ratio evaluations needed
get_log_r <- function(config, idata, pce_vars = NULL, imodel_input, 
                      tmodel_coeffs, proposals, K, N_mixture){

  dim_c <- ncol(tmodel_coeffs) - 2
  
  helper_input <- list(
    S = nrow(tmodel_coeffs),
    N_mixture = N_mixture,
    K = K,
    N_exp = config$n_exp,
    
    y_exp = as.numeric(idata$y_exp),
    
    c_0_mixture = as.numeric(imodel_input$tmodel_coeffs[,"c_0"]),
    c_mixture = imodel_input$tmodel_coeffs[,2:(dim_c + 1), drop = FALSE],
    prob = matrix(rep(1/N_mixture, nrow(tmodel_coeffs)*N_mixture), ncol = N_mixture),
    
    c_0 = as.numeric(tmodel_coeffs[,1]),
    c = tmodel_coeffs[,2:(dim_c+ 1)],
    sigma_sim = as.numeric(tmodel_coeffs[,ncol(tmodel_coeffs)]),
    
    proposal_w = proposals$w_values,
    proposal_sigma = proposals$sigma_values,
    
    ll_weight = 1
  )
  
  # PCE
  if(!is.null(pce_vars)){
    psis_helper_model <- cmdstan_model("stan_code/psis_helper_pce.stan")
    
    helper_input <- append(helper_input,
                           list(M = config$M,
                                d = config$poly_degree,
                                l_poly_coeffs = t(pce_vars$l_poly_coeffs_mat),
                                N_comb = nrow(pce_vars$comb),
                                comb = pce_vars$comb))
  }
  # Logistic
  else{
    psis_helper_model <- cmdstan_model("stan_code/psis_helper_logistic.stan") 
  }
  
  helper_fit <- psis_helper_model$sample(
    data = helper_input,
    seed = config$seed,
    chains = 1,
    parallel_chains = 1,
    iter_sampling = 1,
    iter_warmup = 1,
    refresh = config$refresh_cmdstan,
    adapt_delta = config$adapt_delta,
    fixed_param = TRUE
  )
  
  log_r <- matrix(unlist(helper_fit$draws("log_r", format = "list")), 
                      ncol = nrow(tmodel_coeffs))
  
  # number of log-lik-evaluations
  log_ratio_evals <- (nrow(tmodel_coeffs)+1) * length(proposals$w_values)
  
  return(list("log_r" = log_r,
              "proposals" = proposals,
              "log_ratio_evals" = log_ratio_evals))
}



## Resampling the Posterior distribution using the log ratios
#  Input:
#     - log_weights: matrix[K, S] log ratios from the PSIS
#     - w: matrix[K] draws from the proposal distribution
# Ouptut:
#     - i_posterior: matrix[K, S] draws from the i posterior for each T-draw
resample_psis <- function(log_weights, w){
  S <- ncol(log_weights)
  i_posterior <- matrix(NA, nrow = nrow(log_weights), ncol = ncol(log_weights))
  for(s in 1:S){
    i_posterior[,s] <- sample(w, size = length(w), replace = TRUE,
                              prob = exp(log_weights[,s]))
  }
  return(i_posterior)
}


## Moment Matching
# Input
#   - tmodel_coeffs_MM: tmodel_coeffs of the target values for moment matching
#   - imodel_input_proposal: list: ids, tmodel_coeffs; corresponding to the proposal
#   - imodel: stan-model corresponding to the I-step
#   - t_step: needed for the PCE vars
#   - idata: list with w_exp_1 and y_exp
#   - config
#   - resampling: logical: should resampling be done?
# Output: list
#   - mm: list of moment match objects
#   - log_ratio_evals: number of needed log ratio evaluations
#   - log_prob_evals: number of needed log probability evaluations
#   - pareto_k: pareto k after the moment matching for each imputed dataset
#   - log_weights: log importance weights
#   - proposal: brmsfit object that is used as the proposal
#   - resampled draws: resampled draws using the transformed draws and the log_weights
mm_surrogate <- function(tmodel_coeffs_MM, imodel_input_proposal, 
                         imodel, t_step, idata, config, resampling = FALSE){
  # Log-ratio fun
  log_ratio_fun <- function(draws, fit, tmodel_coeff, proposal_input, pce_vars){
    
    constrained_draws <- rlist_rbind(lapply(apply(
      draws, 1, proposal$constrain_variables), unlist))
    
    dim_c <- ncol(tmodel_coeff) - 2
    helper_input <- list(
      S = nrow(tmodel_coeff),
      N_mixture = 1,
      K = config$K,
      N_exp = config$n_exp,
      
      y_exp = as.numeric(idata$y_exp),
      
      c_0_mixture = as.numeric(proposal_input$tmodel_coeffs[,"c_0"]),
      c_mixture = proposal_input$tmodel_coeffs[,2:(dim_c + 1), drop = FALSE],
      prob = matrix(rep(1/1, nrow(tmodel_coeff)*1), ncol = 1),
      
      c_0 = as.numeric(tmodel_coeff[,1]),
      c = tmodel_coeff[,2:(dim_c+ 1)],
      sigma_sim = as.numeric(tmodel_coeff[,ncol(tmodel_coeff)]),
      
      # Set in new draws
      proposal_w = as.numeric(constrained_draws[,"w_exp"]), # problem: need to be the constrained draws but function takes unconstrained draws
      proposal_sigma = as.numeric(constrained_draws[,"sigma_exp"]),
      
      ll_weight = 1
    )
    
    # PCE
    if(!is.null(pce_vars)){
      psis_helper_model <- cmdstan_model("stan_code/psis_helper_pce.stan")
      helper_input <- append(helper_input,
                             list(M = config$M,
                                  d = config$poly_degree,
                                  l_poly_coeffs = t(pce_vars$l_poly_coeffs_mat),
                                  N_comb = nrow(pce_vars$comb),
                                  comb = pce_vars$comb))
    }
    # Logistic
    else{
      psis_helper_model <- cmdstan_model("stan_code/psis_helper_logistic.stan") 
    }
    
    helper_fit <- psis_helper_model$sample(
      data = helper_input,
      seed = config$seed,
      chains = 1,
      parallel_chains = 1,
      iter_sampling = 1,
      iter_warmup = 1,
      refresh = config$refresh_cmdstan,
      adapt_delta = config$adapt_delta,
      fixed_param = TRUE
    )
    
    log_r <- unlist(helper_fit$draws("log_r", format = "list"))
    
    count_log_ratio <<- count_log_ratio + 2 * nrow(draws)
    
    return(log_r)
  }
  
  # from iwmm:log_prob_draws.CmdStanFit
  log_prob_fun <- function(fit, draws, ...){ # here: draws = unconstrained draws
    count_log_prob <<- count_log_prob + nrow(draws)
    apply(
      draws,
      1,
      fit$log_prob,
      jacobian = TRUE
    )
  }
  
  
  count_log_ratio <<- 0
  count_log_prob <<- 0
  
  # for unconstrained draws: needs to be recompiled in current R Sessions
  x <- try({
    proposal_fit <- fit_imodel(input = imodel_input_proposal, 
                              imodel = imodel, pce_vars = t_step$pce_vars, idata = idata,
                              config = config, ll_weight = 1)
  proposal <- proposal_fit$fits[[1]]
  
  draws <- posterior::as_draws_matrix(proposal)
  udraws <- proposal$unconstrain_draws(format = "draws_matrix")
  })
  if(is(x, "try-error")){
    imodel$compile(force_recompile = TRUE, include_paths = "stan_code") 
    proposal_fit <- fit_imodel(input = imodel_input_proposal, 
                               imodel = imodel, pce_vars = t_step$pce_vars, idata = idata,
                               config = config, ll_weight = 1)
    proposal <- proposal_fit$fits[[1]]
    
    draws <- posterior::as_draws_matrix(proposal)
    udraws <- proposal$unconstrain_draws(format = "draws_matrix")
  }

  mm <- list()
  m <- nrow(tmodel_coeffs_MM)
  
  for(i in 1:m){
    tmodel_coeff <- tmodel_coeffs_MM[i,]
    mm[[i]] <- iwmm:::moment_match.matrix(udraws, log_prob_prop_fun = log_prob_fun,
                                           log_ratio_fun = log_ratio_fun,
                                           fit = proposal,
                                           tmodel_coeff = tmodel_coeff,
                                           proposal_input = proposal_fit$input,
                                           pce_vars = t_step$pce_vars)
  }
  
  mm_draws <- lapply(1:m, function(j) 
    rlist_rbind(lapply(apply(mm[[j]]$draws, 1, function(x) proposal$constrain_variables(x)), unlist)))
  
  for(i in 1:m){
    mm[[i]]$draws <- mm_draws[[i]]
  }
  
  psis_ <- psis(rlist_cbind(lapply(mm, function(x) x$log_weights)))
  log_weights <- psis_$log_weights
  mm_pareto_k <- psis_$diagnostics$pareto_k
  
  # Resampling
  if(resampling){
    draws <- list()
    for(i in 1:m){
      if(mm_pareto_k[i] < 0.7){
        w <- mm_draws[[i]]
        index <- sample(4000, size = length(w), replace = TRUE,
                        prob = exp(psis_$log_weights[,i]))
        draws[[i]] <- as.matrix(fits$fits[[1]])[index,]
      }
    }
  }
  else{
    draws <- NULL
  }
  
  
  
  log_ratio_evals <- count_log_ratio
  log_prob_evals <- count_log_prob
  rm(count_log_ratio, count_log_prob, envir = .GlobalEnv)
  return(list("mm" = mm,
              "log_ratio_evals" = log_ratio_evals,
              "log_prob_evals" = log_prob_evals,
              "pareto_k" = mm_pareto_k,
              "log_weights" = log_weights,
              "proposal" = proposal,
              "resampled_draws" = draws))
  
}

## Helper Functions:
rlist_rbind <- function(list){
  do.call(what = "rbind", args = as.list(list))
}

rlist_cbind <- function(list){
  do.call(what = "cbind", args = as.list(list))
}

