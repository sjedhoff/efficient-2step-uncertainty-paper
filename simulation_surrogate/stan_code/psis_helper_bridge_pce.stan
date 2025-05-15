#include pce_helper.stan

data{
  int<lower=0> S; // Number of t posterior samples
  int<lower=0> N_mixture; // Number of draws used for the mixture proposal
  int<lower=0> K; // number of draws from the proposal for each t posterior draw
  int<lower=0> N_exp; // number of (exp) measured values
  
  // PCE variables
  int<lower=1> M; // dimension of input variable
  int<lower=1> d; // Degree of polynomials
  matrix[d+1, d+1] l_poly_coeffs; // coefficients of legendrepolynomials
  int<lower=1> N_comb; // Number of selected polynomials
  array[N_comb, M] int comb; //polynomial selection
  
  vector[N_exp] y_exp; // Output exp (measured) variables
  
  // T-Posterior draws from the mixture proposal
  vector[N_mixture] c_0_mixture; // values for c_0 from the t-draws selected for the mixture
  matrix[N_mixture, N_comb] c_mixture; // values for c from the t-draws selected for the mixture
  
  // T-Posterior draws of all selected draws
  vector[S] c_0; // values for c_0 from all t-draws
  matrix[S, N_comb] c; // values for c from all t-draws
  vector[S] sigma_sim; // values for sigma_sim all t-draws
  
  // I-Posterior draws for the proposal
  vector[K] proposal_w; // w values from the proposal distribution
  vector[K] proposal_sigma; // sigma values from the proposal distribution
  
  
  vector[N_mixture] log_weights_marginal; // values of the marginal distribution p(y) of the mixture distributions
}

generated quantities{
  matrix[K, S] log_lik_proposal;
  matrix[K, S] log_lik_target;
  matrix[K, S] log_r;

  matrix[K, S] mu_y;
  for(s in 1:S){ // for each t-posterior draw:
  
    // Calculate the mean for the target
    matrix[K,1] w;
    w[:,1] = proposal_w;
    matrix[K, N_comb] pce = get_PCE(w, d, l_poly_coeffs, comb, N_comb);
    mu_y[:,s] = c_0[s] + pce * to_vector(c[s,]);
    
    for(n in 1:K){   // for each w-draw from the proposal
    
          
      // Target Likelihood
      log_lik_target[n, s] =  normal_lpdf(to_vector(y_exp) | mu_y[n,s], proposal_sigma[n]);
      
      
      // Proposal Log-Likelihood    
      vector[1] w_vec;
      w_vec[1] = proposal_w[n]; // select the right w-draw
      
      vector[N_mixture] LSE_input;
      vector[N_mixture] mean_cl;

      real sd = proposal_sigma[n];
      
      // Calculate the mixture using LogSumExp
      for(k in 1:N_mixture){
        real l = 0;
        vector[1] y_exp_pce = c_0_mixture[k] + get_PCE(to_matrix(w_vec,1,1), d, l_poly_coeffs, comb, N_comb) * to_vector(c_mixture[k,]);
        for (i in 1:N_exp) {
          l += normal_lpdf(y_exp[i] | y_exp_pce[1], sd);
        }
        LSE_input[k] = - log_weights_marginal[k] + l;
      }
      log_lik_proposal[n, s] = -log(N_mixture) + log_sum_exp(LSE_input);
      
      // Log-Ratios:
      log_r[n, s] = log_lik_target[n, s] - log_lik_proposal[n, s];

      }
    }
    
}
