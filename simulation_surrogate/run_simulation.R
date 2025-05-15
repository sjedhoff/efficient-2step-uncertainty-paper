################################################################################
##                                                                            ##
##        Simulation: Counting log probability and gradient evaluations       ##
##                for the logistic and the PCE surrogate                      ##
##                                                                            ##
################################################################################


################################################################################
##               First Step: Training the surrogate parameters                ##
################################################################################
source("simulation_surrogate/utils.R")
config <- yaml::read_yaml(file.path("simulation_surrogate/config_sim.yml"), eval.expr=TRUE)

simulator <- function(x1){
  2/(1+exp(-10*(x1-0)))-1
}

## Decide between the logistic and PCE surrogate:
# for the logistic surrogate
t_step <- t_step_logistic(config, simulator)
# for the pce surrogate
t_step <- t_step_pce(config, simulator)



################################################################################
##                         Second Step: Inference                             ##
################################################################################

# for the logistic surrogate
imodel_file <- "simulation_surrogate/stan_code/true_model_logistic_4params_istep_multi_trial.stan"
# for the PCE surrogate
imodel_file <- "simulation_surrogate/stan_code/multi_legendre_pce_istep_mp.stan"

imodel <- cmdstan_model(imodel_file, compile = TRUE)


# simulate the measurement data
idata <- get_idata(config, simulator)

# config gives information about a lot of different hyperparameters, such as
# the number of mixture models etc.. This can be changed like the following:
config$M_mix <- 5

# call the simulation function:
source("simulation_surrogate/simulation_fct_surrogate.R")
result <- simulation_surrogate(t_step = t_step, 
                               imodel = imodel, 
                               idata = idata, 
                               config = config, 
                               method_surrogate = "PCE", # change depending on surrogate type
                               select_method = "random", # choose between maxk, random, llmethod, average
                               method = "PSIS") # change to PSIS or MM (for PSIS+IWMM)
## the results are structured like follows:
# - pareto_k: dataframe with three collumns, where for each event, e.g. PSIS_1 (PSIS in the 
#             first iteration round) for each target_id (t-posterior draws) the corresponding
#             pareto k is saved
# - evals:    dataframe with 6 columns: step (e.g. MCMC_i, PSIS_i, MM_i), log_prob_evals,
#             log_ratio_evals, log_gradient_evals, fits_successful, fits_tried