#include pce_helper.stan

data {
  int<lower=1> N_sim;           // Number of input/output simulation pairs
  int<lower=1> M;               // dimension of input variable
  vector[N_sim] y_sim;          // Output simulation variables
  int<lower=1> d;               // Degree of polynomials
  matrix[N_sim, M] w_sim;       // Input simulation variables

  real sigma_sim_lower;         // Lower bound for sigma_sim
  matrix[d+1, d+1] l_poly_coeffs; // coefficients of legendrepolynomials
  int<lower=1> N_comb; // Number of selected polynomials
  array[N_comb, M] int comb; // polynomial selection
}

transformed data {
  matrix[N_sim, N_comb] W_sim = get_PCE(w_sim, d, l_poly_coeffs, comb, N_comb);
}

parameters {
  real c_0;                                      // coefficient for constant "polynomial"
  vector[N_comb] c;                                   // coefficients of non-constant polynomials
  real<lower=sigma_sim_lower> sigma_sim;         // sigma for simulation values
}


model {
  // Prior model
  c_0 ~ normal(0, 5);
  c ~ normal(0, 5);
  sigma_sim ~ exponential(1);
  
  // Observational model
  y_sim ~ normal_id_glm(W_sim, c_0, c, sigma_sim);
}

generated quantities {
  vector[N_sim] mu_pred_sim;
  mu_pred_sim = W_sim*c + c_0;

  array[N_sim] real y_rep_sim;
  y_rep_sim = normal_rng(mu_pred_sim, sigma_sim);

  vector[N_sim] log_lik_sim;
  for (n in 1:N_sim) log_lik_sim[n] = normal_lpdf(y_sim[n] | mu_pred_sim[n], sigma_sim);
}

