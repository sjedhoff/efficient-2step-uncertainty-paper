#include pce_helper.stan
data {
  int<lower=1> N_exp;          // Number of exp (measured) values
  int<lower=1> N_measures;      // Number of measurements per w_exp
  int<lower=1> M;               // dimension of input variable
  array[N_exp] vector[N_measures] y_exp; // Output exp (measured) variables
  int<lower=1> d;               // Degree of polynomials

  matrix[d+1, d+1] l_poly_coeffs; // coefficients of legendrepolynomials
  int<lower=1> N_comb; // Number of selected polynomials
  array[N_comb, M] int comb; // polynomial selection
  real c_0;                                      // coefficient for constant "polynomial"
  vector[N_comb] c;                                   // coefficients of non-constant polynomials
  
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  real<lower=0> sigma_exp;
  matrix<lower=-1, upper=1>[N_exp, M] w_exp;
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += exponential_lpdf(sigma_exp| 1);
  for (n in 1:N_exp){
    for (m in 1:M){
      lprior += normal_lpdf(w_exp[n, m] | 0, 0.5);
    }
  }
}

model {
  // Observational model
  if (!prior_only) {
    vector[N_exp] y_exp_pce =  c_0 + get_PCE(w_exp, d, l_poly_coeffs, comb, N_comb)*c;
    for (n in 1:N_exp) {
      target +=normal_lpdf(to_vector(y_exp[n, :]) | y_exp_pce[n], sigma_exp);
    }
  }
  // Prior model
  target += lprior;
}

generated quantities{
  vector[N_exp] y_exp_pce =  c_0 + get_PCE(w_exp, d, l_poly_coeffs, comb, N_comb)*c;
  
  vector[N_exp] log_lik;
  for(n in 1:N_exp){
    log_lik[n] = normal_lpdf(to_vector(y_exp[n,:]) | y_exp_pce[n], sigma_exp);
  }
}
