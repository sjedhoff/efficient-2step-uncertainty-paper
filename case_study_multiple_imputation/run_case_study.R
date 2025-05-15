################################################################################
##                                                                            ##
##                    NYC Case Study: Run the simulation                      ##
##                                                                            ##
################################################################################

# run prepare_data.R file first to create the hourly datasets

### Iterate through all hours
for(i in 0:23){
  imp <- readRDS(paste0("case_study_multiple_imputation/data/imp_2016_03_07_",i,".RDS"))
  f <- formula(trip_duration_log ~ osm_distance_scale + osm_duration_scale +
                 pickup_latitude_scale + pickup_longitude_scale + 
                 dropoff_latitude_scale + dropoff_longitude_scale)
  source("simulation_multiple_imputation/utils.R")
  m <- length(imp)
  
  
  ################################################################################
  ##                                  FULL MCMC                                 ##
  ################################################################################
  
  # run MCMC for each imputed dataset separately
  fits <- fit_model(imp = imp, formula = f,
                    marginals = FALSE)
  
  # structure the results and save them
  evals_full_MCMC <- data.frame("step" = paste0("MCMC_", 0),
                                "log_prob_evals" = as.numeric(fits$log_lik_evals[1]),
                                "log_ratio_evals" = 0,
                                "log_gradient_evals" = as.numeric(fits$log_lik_evals[1]),
                                "fits_successful" = m,
                                "fits_tried" = NA)
  total_mcmc_fits <- sum(evals_full_MCMC$fits_successful[grep("MCMC", evals_full_MCMC$step)])
  result_df <- data.frame("dataset" = i,
                          "method" = "MCMC",
                          "total_log_prob_evals" = sum(evals_full_MCMC$log_prob_evals, na.rm = TRUE),
                          "total_log_ratio_evals" = sum(evals_full_MCMC$log_ratio_evals, na.rm = TRUE),
                          "total_log_gradient_evals" = sum(evals_full_MCMC$log_gradient_evals, na.rm = TRUE),
                          "total_MCMC_fits" = total_mcmc_fits,
                          "rounds" = max(as.numeric(gsub(".*_(\\d+)", "\\1", evals_full_MCMC$step)), na.rm = TRUE))
  saveRDS(result_df, file = paste0("case_study_multiple_imputation/results/result_2016_03_07_MCMC_", i, ".rds"))
  
  
  
  ################################################################################
  ##                                  PSIS+MM                                   ##
  ################################################################################
  
  
  algo <- function(imp, N_mix = 1, 
                   select_method = "medoids",
                   method = "MM", 
                   new_mixture = FALSE, 
                   formula = f,
                   ...){
    
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
    
    # Start iteration
    while(length(bad_ids) > N_mix){
      ## Select representative ids
      if(select_method == "maxk" & i > 1){
        selected_id <- bad_ids[which.max(pareto_k_round)]
      }
      if(i == 1){
        selected_id <- sample(bad_ids, N_mix)
      }
      bad_ids <- bad_ids[!(bad_ids %in% selected_id)]
      imp_round <- imp[bad_ids]
      
      ## fit the model to the selected datasets
      fits <- fit_model(imp = imp[selected_id], formula = formula, 
                        marginals = is_mixture, ...)
      
      evals <- rbind(evals, data.frame("step" = paste0("MCMC_", i),
                                       "log_prob_evals" = as.numeric(fits$log_lik_evals[1]),
                                       "log_ratio_evals" = 0,
                                       "log_gradient_evals" = as.numeric(fits$log_lik_evals[1]),
                                       "fits_successful" = length(selected_id),
                                       "fits_tried" = NA))
      
      ## Compute log ratios for PSIS
      log_r <- log_ratios(proposal = fits$fits[[1]], imp = imp_round, is_mixture = is_mixture)
      
      ## Compute PSIS
      psis_ <- psis(log_r$log_r)
      print(psis_)
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
        proposal <- fits$fits[[1]]
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
      remaining_fits <- fit_model(imp = imp[bad_ids], formula = formula, 
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
  
  # Run the algorithm
  result <- algo(imp = imp, N_mix = 1,
                 select_method = "maxk", method = "PSIS", 
                 new_mixture = FALSE, formula = f)
  # save the results
  total_mcmc_fits <- sum(result$evals$fits_successful[grep("MCMC", result$evals$step)])
  result_df <- data.frame("dataset" = i,  
                          "method" = "PSIS",
                          "total_log_prob_evals" = sum(result$evals$log_prob_evals, na.rm = TRUE),
                          "total_log_ratio_evals" = sum(result$evals$log_ratio_evals, na.rm = TRUE),
                          "total_log_gradient_evals" = sum(result$evals$log_gradient_evals, na.rm = TRUE),
                          "total_MCMC_fits" = total_mcmc_fits,
                          "rounds" = max(as.numeric(gsub(".*_(\\d+)", "\\1", result$evals$step)), na.rm = TRUE))
  result[["df"]] <- result_df
  saveRDS(result, file = paste0("case_study_multiple_imputation/results/result_2016_03_07_PSIS_", i, ".rds"))
}