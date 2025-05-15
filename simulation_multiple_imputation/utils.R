################################################################################
##                                                                            ##
##                       Utils for the mutliple Imputation                    ##
##                              Simulation Study                              ##
##                                                                            ##
################################################################################

##### Simulate the data
# Input:
#   - n: number of observations
#   - p: number of features
#   - p_missing: percentage of missing at random
#   - error_variance
# Output:
#   - dataframe (nx(p+1)): variables y and x1, ... xp
library(simstudy)
library(missMethods)
sim_data <- function(n, p, p_missing = 0.1, error_variance = 1){
  
  # make the complete dataset
  v <- rgamma(p+1, 10, 10)
  coef <- sample(x = rnorm(p+1, 0, 1), size = p+1, replace = FALSE)
  C <- matrix(0.3, nrow = p, ncol = p)
  diag(C) <- 1
  data <- genCorData(n = n, mu = numeric(p), sigma = v[1:p], corMatrix = C)
  formula_y <- paste0(coef[1], "+", paste(coef[2:(p+1)], " * V", 1:p, sep = "", collapse = "+"))
  defy <- defDataAdd(varname = "y", formula = formula_y, variance = error_variance, dist = "normal")
  data <- addColumns(defy, data)
  data <- data[,-1]
  data <- as.data.frame(data)[,c((p+1),1:p)]
  
  # fill in missing values
  rows_missing <- sample(1:n, size = ceiling(p_missing * n))
  data[rows_missing, "y"] <- NA
  values_missing <- replicate(length(rows_missing), sample(1:p, size = 0.5 * p, replace = FALSE))
  
  if(!is.matrix(values_missing)){
    values_missing <- as.matrix(t(values_missing))
  }
  
  for(i in 1:length(rows_missing)){
    data[rows_missing[i], values_missing[,i]] <- NA
  }
  
  return(data)
}


##### Imputing Datasets
# Input:
#   - m: number of datasets that should be imputed
#   - data: dataframe, with missing values
# Output:
#   - list with an imputed dataset in each element
library(mice)
impute_data <- function(m, data){
  imp <- mice(data = data, m = m, printFlag = FALSE)
  imp <- complete(imp, "all")
  return(imp)
}


### Pairwise distances between datasets are computed with the Friedman-Rafsky Test
### Input
#   - imp: list with imputed datasets
#   - N_mix: how many datasets should be selected
#   - method: method how to choose the representative datasets: k-medoids or greedy or random
#   - dist: dist matrix, if the distance matrix was calculated before, NULL otherwise
### Output: list
#   - mixture_ids: numeric(N_mix): ids of selected imputed datasets
#   - dist: dist matrix
library(gTests)
library(cluster)
get_imp_representative <- function(imp, N_mix, method = "medoids", dist = NULL){
  
  m <- length(imp)
  if(method == "random"){
    return(list("mixture_ids"= sample(1:m, N_mix),
                "dist" = NULL))
  }
  
  if(is.null(dist)){
    # Friedman-Rafsky-Test
    source("R_multiple_imputation/data_smilarity/helper_functions_gTests.R")
    FR <- function(X1, X2, dist.fun = stats::dist, graph.fun = ade4::mstree, 
                   n.perm = 0, dist.args = NULL, graph.args = NULL, seed = 42) {
      gTestWrapper(X1, X2, dist.fun, graph.fun, n.perm, dist.args, graph.args, 
                   type = "original", seed = seed)
    }
    
    # compute pairwise distances
    dist <- matrix(NA, nrow = m, ncol = m)
    comb <- combn(1:m, 2)
    for(i in 1:ncol(comb)){
      x <- comb[,i]
      dist[x[2], x[1]] <- FR(imp[[x[1]]], imp[[x[2]]])$statistic
    }
    dist <- as.dist(dist)
  }
  
  # k medoids
  if(method == "medoids"){
    k_med <- pam(x = dist, k = N_mix)
    mixture_ids <- k_med$medoids
  }
  # greedy
  else{
    dist_m <- as.matrix(dist)
    # select the one with the maximal mean distance to all other points
    mixture_ids <- which.max(rowMeans(dist_m))
    # get the one with the maximal distance to the first dataset
    mixture_ids[2] <- which.max(dist_m[mixture_ids,])
    # get the next one with the maximal mean distance of the seleted datasets
    for(i in 3:N_mix){
      x <- which.max(colMeans(dist_m[mixture_ids,-mixture_ids]))
      mixture_ids <- c(mixture_ids, as.numeric(names(x)))
    }
  }
  return(list("mixture_ids"= mixture_ids,
              "dist" = dist))
}


##### Fitting models to imputed datasets
# Input:
#   - imp: list of imputed datasets
#   - formula: formula, that needs to be fitted using brms
#   - marginals: logical, should the marginal likelihood be computed using bridge sampling?
# Output: list
#   - fits: list with an brms_fit object for each imputed dataset
#   - log_lik_evals: number of log lik evaluations done in MCMC and bridge sampling, 
#       if num_ll_evals = TRUE
#   - ll_marginals: marginal log likelihood
library(brms)
fit_model <- function(imp, formula, marginals = TRUE, ...){
  m <- length(imp)
  
  n <- nrow(imp[[1]])
  p <- ncol(imp[[1]]) - 1

  
  ## Fitting the models
  fits <- list()
  fits[[1]] <- brm(formula = formula, data = imp[[1]],
                   save_pars = save_pars(all = TRUE),
                   silent = 2, refresh = 0, ...) 
  fits[[1]] <- brms:::update_misc_env(fits[[1]], recompile = FALSE)
  
  ll_marginals <- numeric(m)
  bridges <- list()
  
  if(marginals){
    bridges[[1]] <- bridge_sampler(fits[[1]], silent = TRUE)
    ll_marginals[1] <- bridges[[1]]$logml
  }
  if(m > 1){
    for(i in 2:m){
      fit <- try(update(fits[[1]], newdata = imp[[i]], silent = 2, refresh = 0, ...))
      fits[[i]] <- fit
      fits[[i]] <- brms:::update_misc_env(fits[[i]], recompile = FALSE)
      if(marginals){
        bridges[[i]] <- bridge_sampler(fits[[i]], silent = TRUE)
        ll_marginals[i] <- bridges[[i]]$logml
      }
    }
  }
  
  # calculation of number of log lik evaluations
  log_lik_evals_MCMC <- function(brms_fit) {
    stan_fit <- brms_fit$fit
    sample_params <- rstan::get_sampler_params(stan_fit)
    colnames(sample_params[[1]])
    leapfrog <- lapply(sample_params, function(chain)
      chain[, "n_leapfrog__"])
    log_lik_evals <- sapply(leapfrog, function(f)
      sum(f))
    return(sum(log_lik_evals))
  }
  log_lik_evals_bridge <- function(bridge, fit) {
    bridge$niter * sum(fit$fit@sim$n_save)
  }
  
  log_lik_evals_fit <- unlist(lapply(fits, log_lik_evals_MCMC))
  log_lik_evals_sum <- sum(log_lik_evals_fit)
  if (marginals) {
    log_lik_evals_b <- unlist(lapply(1:m, function(i)
      log_lik_evals_bridge(bridges[[i]], fits[[i]])))
    log_lik_evals_b_sum <- sum(log_lik_evals_b)
    log_lik_evals <- c(MCMC_fits = log_lik_evals_sum, bridge = log_lik_evals_b_sum)
  }
  else{
    log_lik_evals <- c(MCMC_fits = log_lik_evals_sum)
    log_lik_evals_b <- 0
  }
  
  out <- list(
    "fits" = fits,
    "log_lik_evals" = log_lik_evals,
    "ll_marginals" = ll_marginals
  )
  
  return(out)
}



##### Building a mixture
# Input:
#   - fits: list of brms_fits that should be used for the mixture
# Ouput: list
#   - fit: brms_fit object with the mixture draws
#   - draws: draws dataframe of the mixture object with three additional 
#             collumns: draw_id, model_id and post_draw_id
mixture_draws <- function(fits){
  n_sample <- nrow(as.matrix(fits[[1]]))
  N_mix <- length(fits)
  no_draws_model <- sample(1:N_mix, size = n_sample, replace = TRUE)
  draw_ids <- lapply(1:N_mix, function(x) sample(1:n_sample, size = sum(no_draws_model == x), replace = FALSE))
  
  # join into a dataframe
  df_proposal_draw <- data.frame("id" = 1:n_sample,
                                 "model" = no_draws_model,
                                 "post_draw_id" = NA)
  for(i in 1:N_mix){ df_proposal_draw$post_draw_id[df_proposal_draw$model == i] <- draw_ids[[i]]}
  
  # arranging the data frame
  df_proposal_draw <- dplyr::arrange(df_proposal_draw, model, post_draw_id)
  
  # get the posterior values of the selected draws
  values <- lapply(1:N_mix, function(i) as.data.frame(fits[[i]])[df_proposal_draw$post_draw_id[df_proposal_draw$model == i],])
  
  values <- rlist_rbind(values)
  
  df_proposal_draw <- cbind(df_proposal_draw, values)
  
  combine_brmsfit_models <- function(x, draws, ...) {
    draws <- as.matrix(draws)
    ndraws <- nrow(draws)
    ndraws <- ncol(draws)
    stopifnot(all(x$fit@sim$fnames_oi %in% colnames(draws)))
    new_draws <- brms:::named_list(x$fit@sim$fnames_oi, list(numeric(ndraws)))
    for (p in x$fit@sim$fnames_oi) {
      new_draws[[p]] <- draws[, p]
    }
    # create new sim object to overwrite x$fit@sim
    x$fit@sim <- list(
      samples = list(new_draws),
      iter = ndraws,
      thin = 1,
      warmup = 0,
      chains = 1,
      n_save = ndraws,
      warmup2 = 0,
      permutation = list(seq_len(ndraws)),
      pars_oi = x$fit@sim$pars_oi,
      dims_oi = x$fit@sim$dims_oi,
      fnames_oi = x$fit@sim$fnames_oi,
      pars_oi_old = x$fit@sim$pars_oi_old,
      dims_oi_old = x$fit@sim$dims_oi_old,
      fnames_oi_old = x$fit@sim$fnames_oi_old,
      n_flatnames = length(x$fit@sim$fnames_oi)
    )
    x$fit@stan_args <- list(
      list(chain_id = 1, iter = ndraws, thin = 1, warmup = 0)
    )
    x
  }
  mixture <- combine_brmsfit_models(fits[[1]], as.matrix(df_proposal_draw[,-c(1,2,3)]))
  
  return(list("fit" = mixture,
              "draws" = df_proposal_draw))
}



##### Log Ratios
# Input:
#   - proposal: a brms_fit object that represents the proposal distribution
#   - imp: imputed datasets, representing different target distributions
#   - is_mixture: TRUE if the proposal is a mixture, FALSE if its a single proposal
#   - imp_mix: NULL for single proposal, for mixture proposal imputed datasets (of the mixture)
#   - ll_marginals: in case of a mixture proposal, the marginal log likelihoods
# Output: list
#   - log_r: matrix of the log ratios with one column for each imputed dataset
#   - log_lik_evals: number of log_lik evaluations done for computing the log ratios
log_ratios <- function(proposal, imp, is_mixture = FALSE, imp_mix = NULL, ll_marginals = NULL){
  m <- length(imp)
  p <- ncol(proposal$data) - 1
  n <- nrow(proposal$data)
  
  # Target
  diff_rows <- matrix(NA, nrow = m, ncol = n)
  imp_small <- list()
  diff_rows <- rlist_rbind(lapply(1:m, function(i) apply(proposal$data != imp[[i]], 1, any)))
  selected_rows <- apply(diff_rows, 2, any)
  
  # if all proposal = target
  if(sum(selected_rows) == 0){
    log_r <- matrix(1, nrow = nrow(as.matrix(proposal)), ncol = m)
    return(list("log_r" = log_r,
                "log_ratio_evals" = 0))
  }
  for(i in 1:m){
    imp_small[[i]] <- imp[[i]][selected_rows,]
  }
  ll_target <- rlist_cbind(lapply(1:m, function(j) rowSums(log_lik(object = proposal, newdata = imp_small[[j]]))))
  log_ratio_evals_target <- m * nrow(as.matrix(proposal)) * sum(selected_rows)/length(selected_rows)
  
  
  # proposal:
  # single proposal:
  if(!is_mixture){
    ll_proposal <- rowSums(log_lik(object = proposal, newdata = proposal$data[selected_rows,]))
    log_ratio_evals_proposal <- nrow(as.matrix(proposal)) * sum(selected_rows)/length(selected_rows)
  }
  # mixture proposal
  else{
    N_mix <- length(imp_mix)
    if(is.null(ll_marginals)){
      ll_marginals <- rep(0,N_mix)
    }
    
    imp_mix_small <- list()
    for(i in 1:N_mix){
      imp_mix_small[[i]] <- imp_mix[[i]][selected_rows,]
    }
    
    ll_proposal_helper <- rlist_cbind(lapply(1:N_mix, function(j) -ll_marginals[j] + rowSums(log_lik(object = proposal,
                                                                                                     newdata = imp_mix_small[[j]]))))
    ll_proposal <- apply(ll_proposal_helper, 1, matrixStats::logSumExp)
    log_ratio_evals_proposal <- N_mix * nrow(as.matrix(proposal)) * sum(selected_rows)/length(selected_rows)
  }
  
  
  # Log-Ratios:
  log_r <- matrix(NA, nrow = nrow(ll_target), ncol = m)
  for(i in 1:m){
    log_r[,i] <- ll_target[,i] - ll_proposal
    
  }
  
  return(list("log_r" = log_r,
              "log_ratio_evals" = sum(log_ratio_evals_target, log_ratio_evals_proposal)))
}


##### Moment Matching
# Input:
#   - proposal: brmsfit object that is used as the proposal
#   - imp: imputed datasets that are used as the targets
#   - resampling: logical, should the resampling be done?
# Output: list
#   - mm: list of moment match objects
#   - log_ratio_evals: number of needed log ratio evaluations
#   - log_prob_evals: number of needed log probability evaluations
#   - pareto_k: pareto k after the moment matching for each imputed dataset
#   - log_weights: log importance weights
#   - proposal: brmsfit object that is used as the proposal
#   - resampled draws: resampled draws using the transformed draws and the log_weights
### Currently only working for single proposals. For mixture proposals, the following
#   steps could be done: In the mixture_draws-function: brms object with all 
#   N_mix datasets corresponding to the mixture distributions/or containing N_mix
#   corresponding brms-fit/rstan-fit objects. An extra log_prob method
#   needs to be written, which accesses the N_mix datasets/fits and aggregates
#   the log-prob values
mm <- function(proposal, imp, resampling = FALSE){
  # moment match on the brms_fit object using the log_ratio_fun
  log_ratio_fun <- function(draws, fit, imp_data){
    # update the fit with the new draws
    new_fit <- brms:::.update_pars(fit, draws)
    
    # reduce the data to the rows that are different from target vs proposal
    selected_rows <- apply(new_fit$data != imp_data, 1, any)
    imp_data_small <- imp_data[selected_rows,]
    
    # target:
    ll_target <- rowSums(log_lik(object = new_fit, newdata = imp_data_small))
    # proposal
    ll_proposal <- rowSums(log_lik(object = new_fit, newdata = new_fit$data[selected_rows,]))
    # log ratios:
    log_ratios <- ll_target - ll_proposal
    
    count_log_ratio <<- count_log_ratio + 
      nrow(as.matrix(new_fit)) * sum(selected_rows)/length(selected_rows) * 2
    return(log_ratios)
  }
  
  count_log_ratio <<- 0
  count_log_prob <<- 0
  
  log_prob_draws.brmsfit_2 <- function(fit, draws, ...) {
    log_prob_draws.stanfit_2(fit$fit, draws = draws, ...)
  }
  log_prob_draws.stanfit_2 <- function(fit, draws, ...) {
    count_log_prob <<- count_log_prob + nrow(draws)
    apply(
      draws,
      1,
      rstan::log_prob,
      object = fit,
      adjust_transform = TRUE,
      gradient = FALSE
    )
  }
  
  draws <- posterior::as_draws_matrix(proposal)
  udraws <- iwmm:::unconstrain_draws.brmsfit(proposal, draws = draws)
  
  m <- length(imp)
  mm <- lapply(1:m, function(i) iwmm:::moment_match.matrix(udraws,
                                                           log_ratio_fun = log_ratio_fun,
                                                           log_prob_prop_fun = log_prob_draws.brmsfit_2,
                                                           fit = proposal,
                                                           imp_data = imp[[i]]))
  
  mm_draws <- lapply(1:m, function(j) 
    posterior::as_draws_matrix(iwmm:::.update_pars(x = proposal, upars = mm[[j]]$draws)))
  
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
        index <- sample(4000, size = nrow(w), replace = TRUE,
                        prob = exp(psis_$log_weights[,i]))
        draws[[i]] <- w[index,]
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

##### Helper Functions:
rlist_rbind <- function(list){
  do.call(what = "rbind", args = as.list(list))
}

rlist_cbind <- function(list){
  do.call(what = "cbind", args = as.list(list))
}
