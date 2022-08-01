
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed[N_observed];
	int total_observed[N_observed];
  	real<lower=0> rho_hyper_par;
  	real<lower=0> alpha_hyper_par;
}

parameters {
	real<lower=0> rho;
	real<lower=0> alpha;
	real sex0_loc0;
  	vector[N_predict] f_tilde;
}

transformed parameters {
	matrix[N_predict, N_predict] cov;
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict;
	cov = cov_exp_quad(x_predict, alpha, rho) + diag_matrix(rep_vector(1e-10, N_predict));
	L_cov = cholesky_decompose(cov);
	logit_p_predict = sex0_loc0 + L_cov * f_tilde;
}

model {
	rho ~ normal(0, rho_hyper_par);
  	alpha ~ normal(0, alpha_hyper_par);
  	sex0_loc0 ~ normal( 0 , 10 );
  	f_tilde ~ normal(0, 1);
  	y_observed ~ binomial_logit(total_observed, logit_p_predict[observed_idx] );
}

generated quantities {
  vector[N_predict] p_predict = inv_logit(logit_p_predict);  
}
