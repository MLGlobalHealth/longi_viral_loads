
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	vector<lower=0>[N_observed] y_observed;
	vector<lower=0>[N_observed] sd_observed;
	real<lower=0> rho_hyper_par;
	real<lower=0> alpha_hyper_par;
}
			
parameters {
	real<lower=0> rho;
	real<lower=0> alpha;
	real base;
	vector[N_predict] f_tilde;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] mean_predict;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha, rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	mean_predict = base + L_cov * f_tilde;
}
			
model {
	rho ~ normal(0, rho_hyper_par);
	alpha ~ normal(0, alpha_hyper_par);	
	base ~ normal( 0 , 100 );
	f_tilde ~ normal(0, 1);
		
	y_observed ~ normal( mean_predict[observed_idx], sd_observed );
}				
