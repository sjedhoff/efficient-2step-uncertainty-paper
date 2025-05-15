################################################################################
##                                                                            ##
##                  Simulation function for surrogates                        ##
##                                                                            ##
################################################################################
# Input:
#   - t_step: output of t_step function
#   - imodel: compiled cmdstanr model
#   - idata: list with w_exp_1 and y_exp
#   - config: list with configerations
#   - method_surrogate: character, eihter "logistic" or "PCE"
#   - select_method: character: "maxk", "random", "llmethod", "averages"
#   - method: character: "PSIS" or "MM" (PSIS+IWMM)
# Output: List
#   - pareto_k: dataframe with three columns: event, target_id and pareto_k
#   - evals: dataframe counting the number of evaluations in each step:
#       step, log_prob_evals, log_ratio_evals, log_gradient_evals,
#       fits_successful, fits_tried
simulation_surrogate <- function(t_step, imodel, idata, config, method_surrogate = "logistic",
                                 select_method = "random", method = "MM", ...){
  # m different datasets
  m <- nrow(t_step$tmodel_coeffs)
  N_mix <- config$M_mix
  
  # initialize the while loop:
  i <- 1
  is_mixture <- ifelse(N_mix > 1, TRUE, FALSE)
  bad_ids <- 1:m
  length_bad_ids_old <- length(bad_ids)
  
  pareto_k_round <- rep(NA, m)
  evals <- data.frame("step" = character(),
                      "log_prob_evals" = numeric(),
                      "log_ratio_evals" = numeric(),
                      "log_gradient_evals" = numeric(),
                      "fits_successful" = numeric(),
                      "fits_tried" = numeric(),
                      "failed" = logical())
  pareto_k <- data.frame("event" = character(),
                         "target_id" = numeric(),
                         "pareto_k" = numeric())
  
  select_ids <- NULL
  ### Matrix where for each T-Posterior draw (in columns) S draws (in rows) of its I-Posterior will be saved in
  i_posterior <- matrix(NA, nrow = config$K, ncol = config$S, dimnames = list(NULL, 1:config$S))
  
  
  ## for full-MCMC:
  if(method == "MCMC"){
    fits <- fit_imodel(input = list("ids" = 1:m,
                                    "tmodel_coeffs" = t_step$tmodel_coeffs), 
                       imodel = imodel, pce_vars = t_step$pce_vars, idata = idata,
                       config = config, ll_weight = 1)
    evals <- rbind(evals, data.frame("step" = paste0("MCMC_", 0),
                                     "log_prob_evals" = as.numeric(fits$log_lik_evals),
                                     "log_ratio_evals" = 0,
                                     "log_gradient_evals" = as.numeric(fits$log_lik_evals),
                                     "fits_successful" = m,
                                     "fits_tried" = NA))
    return(list("pareto_k" = pareto_k,
                "evals" = evals))
  }
  
  while(length(bad_ids) > N_mix){
    ## Select representative ids
    if(select_method == "maxk" & i > 1){
      if(is_mixture){
        selected_id <-  bad_ids[sort(pareto_k_round, decreasing = TRUE, 
                                     index.return = TRUE, na.last = TRUE)$ix[1:N_mix]]
      }
      else{
        selected_id <- bad_ids[which.max(pareto_k_round)]
      }
      tmodel_coeffs_proposal <- t_step$tmodel_coeffs[selected_id,]
    }
    if(select_method == "random"){
      selected_id <- sample(bad_ids, N_mix)
      tmodel_coeffs_proposal <- t_step$tmodel_coeffs[selected_id,]
    }
    if(select_method == "llmethod" | (i == 1 & select_method == "maxk")){
      selected_t_draws <- select_t_draws(config, t_step$tmodel_coeffs[bad_ids,], idata,
                                         method = method_surrogate)
      selected_id <- selected_t_draws$t_draws_proposal_id
      tmodel_coeffs_proposal <- t_step$tmodel_coeffs[selected_id,]
    }
    if(select_method == "average"){
      selected_id <- NA
      tmodel_coeffs_proposal <- as.matrix(t(apply(t_step$tmodel_coeffs[bad_ids,], 2, mean)))
    }
    
    bad_ids <- bad_ids[!(bad_ids %in% selected_id)]
    tmodel_coeffs_round <- t_step$tmodel_coeffs[bad_ids,]
  
    ## fit the model to the selected datasets -> proposal
    imodel_input_proposal <- list("ids" = selected_id,
                                  "tmodel_coeffs" = tmodel_coeffs_proposal) 
    
    proposal_fit <- fit_imodel(input = imodel_input_proposal, 
               imodel = imodel, pce_vars = t_step$pce_vars, idata = idata,
               config = config, ll_weight = 1)
    if(select_method != "average"){
      i_posterior[,selected_id] <- proposal_fit$w_exp
    }

    N_mixture <- proposal_fit$N_mixture
    # Proposal values
    proposals <- get_proposal_values(w_exp = proposal_fit$w_exp, 
                                     sigma_exp = proposal_fit$sigma_exp)

    evals <- rbind(evals, data.frame("step" = paste0("MCMC_", i),
                                     "log_prob_evals" = as.numeric(proposal_fit$log_lik_evals),
                                     "log_ratio_evals" = 0,
                                     "log_gradient_evals" = as.numeric(proposal_fit$log_lik_evals),
                                     "fits_successful" = length(selected_id),
                                     "fits_tried" = NA))
    
    ## Log-Ratios
    log_r <- get_log_r(config = config, idata = idata, pce_vars = t_step$pce_vars,
                                   imodel_input = proposal_fit$input,
                                   tmodel_coeffs = tmodel_coeffs_round, proposals = proposals,
                                   K = proposal_fit$K, N_mixture = N_mixture)
    
    ## Compute PSIS
    psis_ <- psis(log_r$log_r)
    
    k <- psis_$diagnostics$pareto_k
    
    # set k to 0 if all weights are the same -> target = proposal
    k[apply(psis_$log_weights, 2, function(x) var(x) == 0)] <- 0
    
    evals <- rbind(evals, data.frame("step" = paste0("PSIS_", i),
                                     "log_prob_evals" = 0,
                                     "log_ratio_evals" = as.numeric(log_r$log_ratio_evals),
                                     "log_gradient_evals" = 0,
                                     "fits_successful" = sum(k<= 0.7),
                                     "fits_tried" = length(bad_ids)))
    
    pareto_k <- rbind(pareto_k, data.frame("event" = paste0("PSIS_", i),
                                           "target_id" = bad_ids,
                                           "pareto_k" = k))
    
    pareto_k_round <- k[k > 0.7]
    bad_ids <- bad_ids[k > 0.7]
    if(length(bad_ids) == length_bad_ids_old){
      stop("Average selection not successful")
    }
    length_bad_ids_old <- length(bad_ids)
    
    ## Moment Matching
    if(length(bad_ids) > 0 & method == "MM"){
      if(is_mixture){
        stop("Moment Matching for mixtures is not implemented")
      }
      else{
        tmodel_coeffs_MM <- t_step$tmodel_coeffs[bad_ids,]
        mm_ <- try(mm_surrogate(tmodel_coeffs_MM = t_step$tmodel_coeffs[bad_ids,], 
                            imodel_input_proposal = imodel_input_proposal, 
                            imodel = imodel, t_step = t_step, 
                            idata = idata, config = config, resampling = FALSE))
        if(is(mm_, "try-error")){
          evals <- rbind(evals, data.frame("step" = paste0("MM_", i),
                                           "log_prob_evals" = NA,
                                           "log_ratio_evals" = NA,
                                           "log_gradient_evals" = NA,
                                           "fits_successful" = 0,
                                           "fits_tried" = length(bad_ids)))
          i <- i + 1
          next
        }
        k <- mm_$pareto_k
        pareto_k <- rbind(pareto_k, data.frame("event" = paste0("MM_", i),
                                               "target_id" = bad_ids,
                                               "pareto_k" = k))
        evals <- rbind(evals, data.frame("step" = paste0("MM_", i),
                                         "log_prob_evals" = as.numeric(mm_$log_prob_evals),
                                         "log_ratio_evals" = as.numeric(mm_$log_ratio_evals),
                                         "log_gradient_evals" = 0,
                                         "fits_successful" = sum(k <= 0.7),
                                         "fits_tried" = length(bad_ids)))
        
        bad_ids <- bad_ids[k > 0.7]
        pareto_k_round <- k[k > 0.7]
        
      }

    }
    i <- i + 1
  }
  # fit the remaining with MCMC
  if(length(bad_ids) > 0){
    imodel_input_remaining <- list("ids" = bad_ids,
                                  "tmodel_coeffs" = t_step$tmodel_coeffs[bad_ids,]) 
    
    remaining_fits <- fit_imodel(input = imodel_input_remaining, 
                               imodel = imodel, pce_vars = t_step$pce_vars, idata = idata,
                               config = config, ll_weight = 1)
    
    evals <- rbind(evals, data.frame("step" = paste0("MCMC_", i),
                                     "log_prob_evals" = as.numeric(remaining_fits$log_lik_evals[1]),
                                     "log_ratio_evals" = 0,
                                     "log_gradient_evals" = as.numeric(remaining_fits$log_lik_evals[1]),
                                     "fits_successful" = length(bad_ids),
                                     "fits_tried" = NA))
  }
  
  return(list("pareto_k" = pareto_k,
              "evals" = evals))
  
  
}
