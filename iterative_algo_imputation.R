################################################################################
#                                                                              #
#             Iterative Algorithm for linear models with missing data          #
#                                                                              #
################################################################################

library(mice)
library(brms)
### Iterative Algorithm
## Input:
#   - data: dataframe of your choice
#   - formula: formula for modeling variables within your dataset
#   - m: number of imputated dataset constructed with MICE
#   - ...: additional inputs for the brm function
## Output:
#   - draws: each list element contains a draws matrix corresponding to the imputed datasets
#   - imp: list of the m imputed datasets
iterative_algo_imputation <- function(data, formula, m = 20, ...){
  
  ################################################################################
  ##                            Helper Functions                                ##
  ################################################################################
  
  
  ##### Log Ratios
  # Input:
  #   - proposal: a brms_fit object that represents the proposal distribution
  #   - imp: imputed datasets, representing different target distributions
  # Output: list
  #   - log_r: matrix of the log ratios with one column for each imputed dataset
  #   - log_lik_evals: number of log_lik evaluations done for computing the log ratios
  log_ratios <- function(proposal, imp){
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
    ll_proposal <- rowSums(log_lik(object = proposal, newdata = proposal$data[selected_rows,]))
    log_ratio_evals_proposal <- nrow(as.matrix(proposal)) * sum(selected_rows)/length(selected_rows)
    
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
  rlist_rbind <- function(list){
    do.call(what = "rbind", args = as.list(list))
  }
  
  rlist_cbind <- function(list){
    do.call(what = "cbind", args = as.list(list))
  }
  
  ##############################################################################
  ##                                 Algorithm                                ##
  ##############################################################################
  
  # impute the dataset m times using the mice package
  imp_all <- mice(data = data, m = m, printFlag = FALSE)
  imp_all <- complete(imp_all, "all")
  
  # Filter variables which are needed in the formula
  imp <- lapply(imp_all, function(x) x[,all.vars(formula)])
  
  bad_ids <- 1:m
  select_ids <- NULL
  
  draws <- list()
  selected_id <- sample(bad_ids,1)
  
  while(length(bad_ids) > 1){
    bad_ids <- bad_ids[!(bad_ids %in% selected_id)]
    imp_round <- imp[bad_ids]
    
    # fit the model using MCMC on the proposal
    proposal <- brm(formula = formula, data = imp[[selected_id]],
                    save_pars = save_pars(all = TRUE),
                    silent = 2, refresh = 0, ...)
    draws[[selected_id]] <- as_draws_matrix(proposal)
    
    # calculate importance ratios and run PSIS
    log_r <- log_ratios(proposal = proposal, imp = imp_round)
    psis_ <- psis(log_r$log_r)
    k <- psis_$diagnostics$pareto_k
    
    # Importance weighted Resampling
    for(i in bad_ids){
      if(k[which(bad_ids %in% i)] < 0.7){
        w <- draws[[selected_id]]
        index <- sample(4000, size = nrow(w), replace = TRUE,
                        prob = exp(psis_$log_weights[,which(bad_ids %in% i)]))
        draws[[i]] <- w[index,]
      }
    }
    
    pareto_k_round <- k[k > 0.7]
    bad_ids <- bad_ids[k > 0.7]
    
    if(length(bad_ids) > 0){
      # IWMM
      imp_MM <- imp[bad_ids]
      mm_ <- mm(proposal = proposal, imp = imp_MM, resampling = TRUE)
      
      k <- mm_$pareto_k
      draws[bad_ids[k < 0.7]] <- mm_$resampled_draws
      
      bad_ids <- bad_ids[k > 0.7]
      pareto_k_round <- k[k > 0.7]
    }
    
    # chose a representative for the proposal
    selected_id <- bad_ids[which.max(pareto_k_round)]
  }
  return(list("draws" = draws,
              "imp" = imp))
}





