################################################################################
##                                                                            ##
##          Simulation function for the multiple imputation case              ##
##                                                                            ##
################################################################################

# Input:
#   - imp: list with imputed datasets in each list entry.
#   - N_mix: number of distributions used for the mixture proposal. 
#       N_mix = 1 means no mixture proposal.
#   - select_method: character: "random", "maxk" or "medoids"
#   - method: character: "PSIS" or "MM" (PSIS+IWMM, only for N_mix = 1)
#   - new_mixture: logical: TRUE or FALSE. If N_mix > 1: should a new
#       mixture be build in each iteration or should the mixture be extended
# Output: List
#   - pareto_k: dataframe with three columns: event, target_id and pareto_k
#   - evals: dataframe counting the number of evaluations in each step:
#       step, log_prob_evals, log_ratio_evals, log_gradient_evals,
#       fits_successful, fits_tried
simulation_imputation <- function(imp, N_mix = 1,
                       select_method = "random", method = "PSIS", 
                       new_mixture = FALSE, ...){
  
  # m different datasets
  m <- length(imp)
  
  # initialize the while loop:
  i <- 1
  is_mixture <- ifelse(N_mix > 1, TRUE, FALSE)
  bad_ids <- 1:m
  
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
  
  ## for full-MCMC:
  if(method == "MCMC"){
    fits <- fit_model(imp = imp, formula = as.formula("y ~ .", env = new.env()), 
                      marginals = TRUE, ...)
    evals <- rbind(evals, data.frame("step" = paste0("MCMC_", 0),
                                     "log_prob_evals" = as.numeric(fits$log_lik_evals[1]),
                                     "log_ratio_evals" = 0,
                                     "log_gradient_evals" = as.numeric(fits$log_lik_evals[1]),
                                     "fits_successful" = m,
                                     "fits_tried" = NA))
    return(list("pareto_k" = pareto_k,
                "evals" = evals))
  }
  
  
  # Start iteration
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
    }
    if(select_method == "random"){
      selected_id <- sample(bad_ids, N_mix)
    }
    if(select_method == "medoids" | (i == 1 & select_method == "maxk")){
      # calculate the complete distance matrix for the first run
      if(i == 1){
        total_dist <- get_imp_representative(imp = imp,
                                             N_mix = N_mix,
                                             method = "medoids")$dist
      }
      
      part_dist <- as.dist(as.matrix(total_dist)[bad_ids, bad_ids])
      selected_id <- get_imp_representative(imp = imp[bad_ids],
                                            N_mix = N_mix,
                                            method = "medoids",
                                            dist = part_dist)$mixture_ids
      selected_id <- as.numeric(selected_id)
      
    }
    
    bad_ids <- bad_ids[!(bad_ids %in% selected_id)]
    imp_round <- imp[bad_ids]
    
    ## fit the model to the selected datasets
    fits <- fit_model(imp = imp[selected_id], formula = as.formula("y ~ .", env = new.env()),
                      marginals = is_mixture, ...)
    
    evals <- rbind(evals, data.frame("step" = paste0("MCMC_", i),
                                     "log_prob_evals" = as.numeric(fits$log_lik_evals[1]),
                                     "log_ratio_evals" = 0,
                                     "log_gradient_evals" = as.numeric(fits$log_lik_evals[1]),
                                     "fits_successful" = length(selected_id),
                                     "fits_tried" = NA))
    if(is_mixture){
      evals <- rbind(evals, data.frame("step" = paste0("bridge_",i),
                                       "log_prob_evals" = as.numeric(ifelse(is.na(fits$log_lik_evals[2]), 0, fits$log_lik_evals[2])),
                                       "log_ratio_evals" = 0,
                                       "log_gradient_evals" = 0,
                                       "fits_successful" = 0,
                                       "fits_tried" = NA))
    }
    
    ## Compute log ratios for PSIS
    # for mixtures
    if(is_mixture){
      mixture <- mixture_draws(fits = fits$fits)
      log_r <- log_ratios(proposal = mixture$fit, imp = imp_round, is_mixture = is_mixture,
                          imp_mix = imp[selected_id], ll_marginals = fits$ll_marginals)
    }
    # single
    else{
      log_r <- log_ratios(proposal = fits$fits[[1]], imp = imp_round, is_mixture = is_mixture)
    }
    
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
    
    ## Moment Matching
    if(length(bad_ids) > 0 & method == "MM"){
      if(is_mixture){
        stop("Moment Matching for mixtures is not implemented")
      }
      else{
        proposal <- fits$fits[[1]]
      }
      imp_MM <- imp[bad_ids]
      proposal <- brms:::update_misc_env(proposal, recompile = FALSE)
      mm_ <- try(mm(proposal = proposal, imp = imp_MM))
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
    i <- i + 1
  } # end while
  
  # fit the remaining with MCMC
  if(length(bad_ids) > 0){
    remaining_fits <- fit_model(imp = imp[bad_ids], formula = as.formula("y ~ .", env = new.env()), 
                                marginals = FALSE, ...)
    
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

