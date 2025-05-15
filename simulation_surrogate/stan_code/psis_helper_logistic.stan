functions {
  real logistic_fct_mean_vec(real c_0, vector c, real w) {
    return c_0 / (1 + exp(-c[1] * (w - c[2]))) + c[3];
  }
  
}


data{
  int<lower=0> S; // Number of t posterior samples
  int<lower=0> N_mixture; // Number of draws used for the mixture proposal
  int<lower=0> K; // number of draws from the proposal for each t posterior draw
  int<lower=0> N_exp; // number of (exp) measured values
  
  vector[N_exp] y_exp; // Output exp (measured) variables
  
  // T-Posterior draws from the mixture proposal
  vector[N_mixture] c_0_mixture; // values for c_0 from the t-draws selected for the mixture
  matrix[N_mixture, 3] c_mixture; // values for c from the t-draws selected for the mixture
  matrix[S, N_mixture] prob; // probability of each t posterior draw corresponding to the clusters
  
  // T-Posterior draws of all selected draws
  vector[S] c_0; // values for c_0 from all t-draws
  matrix[S, 3] c; // values for c from all t-draws
  vector[S] sigma_sim; // values for sigma_sim all t-draws
  
  // I-Posterior draws for the proposal
  vector[K] proposal_w; // w values from the proposal distribution
  vector<lower=0>[K] proposal_sigma; // sigma values from the proposal distribution
  
  
  real ll_weight; // log likelihood weight
}

generated quantities{
  matrix[K, S] log_lik_proposal;
  matrix[K, S] log_lik_target;
  matrix[K, S] log_r;
  
  matrix[K, S] mu_y;
  for(s in 1:S){ // for each t posterior draw:
  
    // Calculate the mean
    mu_y[:,s] = c_0[s] / (1 + exp(-c[s,1] * (proposal_w - c[s,2]))) + c[s,3];

    for(n in 1:K){   // for each w-draw from the proposal
      // Target Likelihood
      log_lik_target[n, s] = ll_weight * normal_lpdf(to_vector(y_exp) | mu_y[n,s], proposal_sigma[n]);
      
      // Proposal Likelihood
      real w = proposal_w[n];
      vector[N_mixture] LSE_input;
      vector[N_mixture] mean_cl;

      real sd = proposal_sigma[n];
      
      // Calculate the mixture using LogSumExp
      for(k in 1:N_mixture){
        real l = 0;
        for(i in 1:N_exp){
          l += ll_weight * normal_lpdf(y_exp[i] | logistic_fct_mean_vec(c_0_mixture[k], to_vector(c_mixture[k]), w), sd);
        }
        LSE_input[k] = log(prob[s,k]) + l;
      }
    
      log_lik_proposal[n, s] = log_sum_exp(LSE_input);
      
      // Log-Ratios:
      log_r[n, s] = log_lik_target[n, s] - log_lik_proposal[n, s];
    }

  }
}
