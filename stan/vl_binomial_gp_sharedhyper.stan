// Changes:
// - hyperprior for lengthscale set to inverse gamma
// - specified by informind lower and upper bound st P(outside bounds)==.02

functions {
  vector tail_delta(vector y, vector theta, array[] real x_r, array[] int x_i) {
    vector[2] deltas;
    deltas[1] = inv_gamma_cdf(theta[1] | exp(y[1]), exp(y[2])) - 0.01;
    deltas[2] = 1 - inv_gamma_cdf(theta[2] | exp(y[1]), exp(y[2])) - 0.01;
    return deltas;
  }
  
  vector make_logits_from_gp(int N_predict, array[] real x_predict,
                             real alpha, real rho, real baseline,
                             vector f_tilde) {
    matrix[N_predict, N_predict] L_cov;
    vector[N_predict] logit_p_predict;
    
    L_cov = cholesky_decompose(gp_exp_quad_cov(x_predict, alpha, rho)
                               + diag_matrix(rep_vector(1e-10, N_predict)));
    logit_p_predict = baseline + L_cov * f_tilde;
    
    return logit_p_predict;
  }
}
data {
  //  number of age-groups
  int<lower=1> N_predict;
  array[N_predict] real x_predict;
  int<lower=1> N_observed;
  array[N_observed] int<lower=1, upper=N_predict> observed_idx;
  
  // numerators
  array[N_observed] int y_observed_00;
  array[N_observed] int y_observed_10;
  array[N_observed] int y_observed_01;
  array[N_observed] int y_observed_11;
  // denominators
  array[N_observed] int total_observed_00;
  array[N_observed] int total_observed_10;
  array[N_observed] int total_observed_01;
  array[N_observed] int total_observed_11;
  
  // First-time participants
  // numerators 
  array[N_observed] int y_observed_00_ftp;
  array[N_observed] int y_observed_10_ftp;
  array[N_observed] int y_observed_01_ftp;
  array[N_observed] int y_observed_11_ftp;
  // denominators
  array[N_observed] int total_observed_00_ftp;
  array[N_observed] int total_observed_10_ftp;
  array[N_observed] int total_observed_01_ftp;
  array[N_observed] int total_observed_11_ftp;
  
  // these GP hyperparameters are shared between ftp and non-ftp across each 
  // (female-male) x (fishing-inland) groups
  real<lower=0> rho_hyper_lower_bound;
  real<lower=0> rho_hyper_upper_bound;
  real<lower=0> alpha_hyper_par_00;
  real<lower=0> alpha_hyper_par_10;
  real<lower=0> alpha_hyper_par_01;
  real<lower=0> alpha_hyper_par_11;
}
transformed data {
  // https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html
  vector[2] log_guess = [log(10), log(20)]';
  vector[2] theta = [rho_hyper_lower_bound, rho_hyper_upper_bound]';
  vector[2] y;
  real<lower=0> rho_hyper_par_shape;
  real<lower=0> rho_hyper_par_scale;
  array[0] real x_r;
  array[0] int x_i;
  
  y = algebra_solver(tail_delta, log_guess, theta, x_r, x_i);
  rho_hyper_par_shape = exp(y[1]);
  rho_hyper_par_scale = exp(y[2]);
  
  print("prior shape = ", exp(y[1]));
  print("prior scale = ", exp(y[2]));
}
parameters {
  real<lower=0> rho_00;
  real<lower=0> rho_10;
  real<lower=0> rho_01;
  real<lower=0> rho_11;
  real<lower=0> alpha_00;
  real<lower=0> alpha_10;
  real<lower=0> alpha_01;
  real<lower=0> alpha_11;
  real sex0_loc0;
  real sex1_loc0;
  real sex0_loc1;
  real sex1_loc1;
  real sex0_loc0_ftp;
  real sex1_loc0_ftp;
  real sex0_loc1_ftp;
  real sex1_loc1_ftp;
  vector[N_predict] f_tilde_00;
  vector[N_predict] f_tilde_10;
  vector[N_predict] f_tilde_01;
  vector[N_predict] f_tilde_11;
  vector[N_predict] f_tilde_00_ftp;
  vector[N_predict] f_tilde_10_ftp;
  vector[N_predict] f_tilde_01_ftp;
  vector[N_predict] f_tilde_11_ftp;
}
transformed parameters {
  // matrix[N_predict, N_predict] L_cov;
  vector[N_predict] logit_p_predict_00;
  vector[N_predict] logit_p_predict_10;
  vector[N_predict] logit_p_predict_01;
  vector[N_predict] logit_p_predict_11;
  vector[N_predict] logit_p_predict_00_ftp;
  vector[N_predict] logit_p_predict_10_ftp;
  vector[N_predict] logit_p_predict_01_ftp;
  vector[N_predict] logit_p_predict_11_ftp;
  
  logit_p_predict_00 = make_logits_from_gp(N_predict, x_predict, alpha_00,
                                           rho_00, sex0_loc0, f_tilde_00);
  logit_p_predict_01 = make_logits_from_gp(N_predict, x_predict, alpha_01,
                                           rho_01, sex0_loc1, f_tilde_01);
  logit_p_predict_10 = make_logits_from_gp(N_predict, x_predict, alpha_10,
                                           rho_10, sex1_loc0, f_tilde_10);
  logit_p_predict_11 = make_logits_from_gp(N_predict, x_predict, alpha_11,
                                           rho_11, sex1_loc1, f_tilde_11);
  
  logit_p_predict_00_ftp = make_logits_from_gp(N_predict, x_predict,
                                               alpha_00, rho_00,
                                               sex0_loc0_ftp, f_tilde_00_ftp);
  logit_p_predict_01_ftp = make_logits_from_gp(N_predict, x_predict,
                                               alpha_01, rho_01,
                                               sex0_loc1_ftp, f_tilde_01_ftp);
  logit_p_predict_10_ftp = make_logits_from_gp(N_predict, x_predict,
                                               alpha_10, rho_10,
                                               sex1_loc0_ftp, f_tilde_10_ftp);
  logit_p_predict_11_ftp = make_logits_from_gp(N_predict, x_predict,
                                               alpha_11, rho_11,
                                               sex1_loc1_ftp, f_tilde_11_ftp);
}
model {
  rho_00 ~ inv_gamma(rho_hyper_par_shape, rho_hyper_par_scale);
  rho_10 ~ inv_gamma(rho_hyper_par_shape, rho_hyper_par_scale);
  rho_01 ~ inv_gamma(rho_hyper_par_shape, rho_hyper_par_scale);
  rho_11 ~ inv_gamma(rho_hyper_par_shape, rho_hyper_par_scale);
  alpha_00 ~ normal(0, alpha_hyper_par_00);
  alpha_10 ~ normal(0, alpha_hyper_par_10);
  alpha_01 ~ normal(0, alpha_hyper_par_01);
  alpha_11 ~ normal(0, alpha_hyper_par_11);
  
  sex0_loc0 ~ normal(0, 10);
  sex0_loc1 ~ normal(0, 10);
  sex1_loc0 ~ normal(0, 10);
  sex1_loc1 ~ normal(0, 10);
  f_tilde_00 ~ normal(0, 1);
  f_tilde_01 ~ normal(0, 1);
  f_tilde_10 ~ normal(0, 1);
  f_tilde_11 ~ normal(0, 1);
  sex0_loc0_ftp ~ normal(0, 10);
  sex0_loc1_ftp ~ normal(0, 10);
  sex1_loc0_ftp ~ normal(0, 10);
  sex1_loc1_ftp ~ normal(0, 10);
  f_tilde_00_ftp ~ normal(0, 1);
  f_tilde_01_ftp ~ normal(0, 1);
  f_tilde_10_ftp ~ normal(0, 1);
  f_tilde_11_ftp ~ normal(0, 1);
  y_observed_00 ~ binomial_logit(total_observed_00,
                                 logit_p_predict_00[observed_idx]);
  y_observed_01 ~ binomial_logit(total_observed_01,
                                 logit_p_predict_01[observed_idx]);
  y_observed_10 ~ binomial_logit(total_observed_10,
                                 logit_p_predict_10[observed_idx]);
  y_observed_11 ~ binomial_logit(total_observed_11,
                                 logit_p_predict_11[observed_idx]);
  
  y_observed_00_ftp ~ binomial_logit(total_observed_00_ftp,
                                     logit_p_predict_00_ftp[observed_idx]);
  y_observed_01_ftp ~ binomial_logit(total_observed_01_ftp,
                                     logit_p_predict_01_ftp[observed_idx]);
  y_observed_10_ftp ~ binomial_logit(total_observed_10_ftp,
                                     logit_p_predict_10_ftp[observed_idx]);
  y_observed_11_ftp ~ binomial_logit(total_observed_11_ftp,
                                     logit_p_predict_11_ftp[observed_idx]);
}
generated quantities {
  // transformed to probability space
  vector[N_predict] p_predict_00;
  vector[N_predict] p_predict_01;
  vector[N_predict] p_predict_10;
  vector[N_predict] p_predict_11;
  
  vector[N_predict] p_predict_00_ftp;
  vector[N_predict] p_predict_01_ftp;
  vector[N_predict] p_predict_10_ftp;
  vector[N_predict] p_predict_11_ftp;
  // output hyperparameters
  real rho_hyper_par_shape2;
  real rho_hyper_par_scale2;
  // // log-likelihood
  array[N_observed] real log_lik_00;
  array[N_observed] real log_lik_01;
  array[N_observed] real log_lik_10;
  array[N_observed] real log_lik_11;
  // // samples from PPD
  // array[N_predict] int y_observed_00_rep;
  // array[N_predict] int y_observed_01_rep;
  // array[N_predict] int y_observed_10_rep;
  // array[N_predict] int y_observed_11_rep;
  // 
  // // compute pointwise log likelihood and posterior predictive at every iteration
  // 
  for (i in 1 : N_observed) {
    log_lik_00[i] = binomial_logit_lpmf(y_observed_00[i] | total_observed_00[i], logit_p_predict_00[i]);
    log_lik_01[i] = binomial_logit_lpmf(y_observed_01[i] | total_observed_01[i], logit_p_predict_01[i]);
    log_lik_10[i] = binomial_logit_lpmf(y_observed_10[i] | total_observed_10[i], logit_p_predict_10[i]);
    log_lik_11[i] = binomial_logit_lpmf(y_observed_11[i] | total_observed_11[i], logit_p_predict_11[i]);
  //   
  //   // y_observed_00_rep[i] ~ binomial_logit(total_observed_00[i], logit_p_predict_00[i]);
  //   // y_observed_01_rep[i] ~ binomial_logit(total_observed_01[i], logit_p_predict_01[i]);
  //   // y_observed_10_rep[i] ~ binomial_logit(total_observed_10[i], logit_p_predict_10[i]);
  //   // y_observed_11_rep[i] ~ binomial_logit(total_observed_11[i], logit_p_predict_11[i]);
  }
  
  p_predict_00 = inv_logit(logit_p_predict_00);
  p_predict_01 = inv_logit(logit_p_predict_01);
  p_predict_10 = inv_logit(logit_p_predict_10);
  p_predict_11 = inv_logit(logit_p_predict_11);
  p_predict_00_ftp = inv_logit(logit_p_predict_00_ftp);
  p_predict_01_ftp = inv_logit(logit_p_predict_01_ftp);
  p_predict_10_ftp = inv_logit(logit_p_predict_10_ftp);
  p_predict_11_ftp = inv_logit(logit_p_predict_11_ftp);
  rho_hyper_par_shape2 = rho_hyper_par_shape;
  rho_hyper_par_scale2 = rho_hyper_par_scale;
}

