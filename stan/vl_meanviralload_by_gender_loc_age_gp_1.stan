
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	int zero_observed[N_observed];
	int total_observed[N_observed];
	vector<lower=0>[N_observed] mean_observed;
	vector<lower=0>[N_observed] meansd_observed;
	real<lower=0> zero_rho_hyper_par;
	real<lower=0> zero_alpha_hyper_par;
	real<lower=0> mean_rho_hyper_par;
	real<lower=0> mean_alpha_hyper_par;
}
			
parameters {
	real<lower=0> zero_rho;
	real<lower=0> zero_alpha;
	real zero_base;
	vector[N_predict] zero_f_tilde;
	real<lower=0> mean_rho;
	real<lower=0> mean_alpha;			
	real mean_base;
	vector[N_predict] mean_f_tilde;
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] logit_zero_p_predict;
	vector[N_predict] zero_p_predict;
	vector[N_predict] mean_predict;
	// GP on zeros
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, zero_alpha, zero_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_zero_p_predict = zero_base + L_cov * zero_f_tilde;
	zero_p_predict = inv_logit(logit_zero_p_predict);
	// GP on means
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, mean_alpha, mean_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	mean_predict = mean_base + L_cov * mean_f_tilde;

}
			
model {
	zero_rho ~ normal(0, zero_rho_hyper_par);
	zero_alpha ~ normal(0, zero_alpha_hyper_par);
	zero_base ~ normal( 0 , 100 );
	zero_f_tilde ~ normal(0, 1);
	
	mean_rho ~ normal(0, mean_rho_hyper_par);
	mean_alpha ~ normal(0, mean_alpha_hyper_par);
	mean_base ~ normal( 0 , 100 );
	mean_f_tilde ~ normal(0, 1);
	
	zero_observed ~ binomial_logit( total_observed, logit_zero_p_predict[observed_idx] );
	mean_observed ~ normal( (1-zero_p_predict[observed_idx]) .* mean_predict[observed_idx], meansd_observed + rep_vector(1e-10, N_observed));
}
			
generated quantities {
	vector[N_predict] zero_inflated_mean_predict = (1-zero_p_predict) .* mean_predict;  
}			

