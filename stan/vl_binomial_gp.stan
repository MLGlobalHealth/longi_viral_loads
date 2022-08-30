// Changes:
// - hyperprior for lengthscale set to inverse gamma
// - specified by informind lower and upper bound st P(outside bounds)==.02

functions{
  vector tail_delta(vector y, vector theta, real[] x_r, int[] x_i) {
    vector[2] deltas;
    deltas[1] = inv_gamma_cdf(theta[1], exp(y[1]), exp(y[2])) - 0.01;
    deltas[2] = 1 - inv_gamma_cdf(theta[2], exp(y[1]), exp(y[2])) - 0.01;
    return deltas;
    }
}

data{	
	int<lower=1> N_predict;
  real x_predict[N_predict];
  int<lower=1> N_observed;
  int<lower=1, upper=N_predict> observed_idx[N_observed];
  int y_observed_00[N_observed];
	int y_observed_10[N_observed];
	int y_observed_01[N_observed];	
	int y_observed_11[N_observed];	
	int total_observed_00[N_observed];
	int total_observed_10[N_observed];
	int total_observed_01[N_observed];
	int total_observed_11[N_observed];
  real<lower=0> rho_hyper_lower_bound;
  real<lower=0> rho_hyper_upper_bound;
  real<lower=0> alpha_hyper_par_00;
	real<lower=0> alpha_hyper_par_10;
	real<lower=0> alpha_hyper_par_01;
	real<lower=0> alpha_hyper_par_11;
}

transformed data{
  // https://betanalpha.github.io/assets/case_studies/gp_part3/part3.html
  vector[2] log_guess = [log(10), log(20)]';
  vector[2] theta = [rho_hyper_lower_bound, rho_hyper_upper_bound]';
  vector[2] y;
  real<lower=0> rho_hyper_par_shape;
  real<lower=0> rho_hyper_par_scale;
  real x_r[0];
  int x_i[0];

  y = algebra_solver(tail_delta, log_guess, theta, x_r, x_i);
  rho_hyper_par_shape = exp( y[1] );
  rho_hyper_par_scale = exp( y[2] );

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
  vector[N_predict] f_tilde_00;
	vector[N_predict] f_tilde_10;
	vector[N_predict] f_tilde_01;
	vector[N_predict] f_tilde_11;
}

transformed parameters {
  matrix[N_predict, N_predict] L_cov;
  vector[N_predict] logit_p_predict_00;
	vector[N_predict] logit_p_predict_10;
	vector[N_predict] logit_p_predict_01;
	vector[N_predict] logit_p_predict_11;
	// GP for 00 and 01 (women)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_00, rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_00 = sex0_loc0 + L_cov * f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_01, rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_01 = sex0_loc1 + L_cov * f_tilde_01;
	// GP for 10 and 10 (men)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_10, rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_10 = sex1_loc0 + L_cov * f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_11, rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_11 = sex1_loc1 + L_cov * f_tilde_11;
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
  sex0_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
  f_tilde_00 ~ normal(0, 1);
	f_tilde_01 ~ normal(0, 1);
	f_tilde_10 ~ normal(0, 1);
	f_tilde_11 ~ normal(0, 1);
  y_observed_00 ~ binomial_logit(total_observed_00, logit_p_predict_00[observed_idx] );
	y_observed_01 ~ binomial_logit(total_observed_01, logit_p_predict_01[observed_idx] );
	y_observed_10 ~ binomial_logit(total_observed_10, logit_p_predict_10[observed_idx] );
	y_observed_11 ~ binomial_logit(total_observed_11, logit_p_predict_11[observed_idx] );
}

generated quantities {
  vector[N_predict] p_predict_00;
	vector[N_predict] p_predict_01;
	vector[N_predict] p_predict_10;
	vector[N_predict] p_predict_11;
  real rho_hyper_par_shape2;
  real rho_hyper_par_scale2;

	p_predict_00 = inv_logit(logit_p_predict_00);
	p_predict_01 = inv_logit(logit_p_predict_01);  
	p_predict_10 = inv_logit(logit_p_predict_10);
	p_predict_11 = inv_logit(logit_p_predict_11);
  rho_hyper_par_shape2=rho_hyper_par_shape;
  rho_hyper_par_scale2=rho_hyper_par_scale; 
}


