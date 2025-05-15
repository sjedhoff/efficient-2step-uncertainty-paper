################################################################################
##                                                                            ##
##        Simulation: Counting log probability and gradient evaluations       ##
##                   for the Multiple Imputation Setup                        ##
##                                                                            ##
################################################################################

################################################################################
##            First Step: Generate Dataset and impute it with MICE            ##
################################################################################

source("simulation_multiple_imputation/utils.R")
# simulate a dataset with size n, predictors p, and percentage of rows for missing values
data <- sim_data(n = 100, p = 20, p_missing = 0.15)
imp <- impute_data(m = 100, data = data)

################################################################################
##                         Second Step: Inference                             ##
################################################################################

source("simulation_multiple_imputation/simulation_fct.R")
result <- simulation_imputation(imp = imp, N_mix = 5,
                                select_method = "random", # select one of maxk, random, llmethod, average
                                method = "PSIS", # "PSIS" or "MM" (PSIS+IWMM)
                                new_mixture = FALSE)
## the results are structured like follows:
# - pareto_k: dataframe with three collumns, where for each event, e.g. PSIS_1 (PSIS in the 
#             first iteration round) for each target_id (imputed dataset) the corresponding
#             pareto k is saved
# - evals:    dataframe with 6 columns: step (e.g. MCMC_i, PSIS_i, MM_i), log_prob_evals,
#             log_ratio_evals, log_gradient_evals, fits_successful, fits_tried