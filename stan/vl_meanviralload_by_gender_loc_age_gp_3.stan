
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	vector<lower=0>[N_observed] y_observed;
	vector<lower=0>[N_observed] s2_observed; // with denominator (n-1)
	vector<lower=1>[N_observed] n_observed;
	real<lower=0> m_rho_hyper_par_alpha;
	real<lower=0> m_rho_hyper_par_beta;
	real<lower=0> m_alpha_hyper_par;
	real<lower=0> s_rho_hyper_par_alpha;
	real<lower=0> s_rho_hyper_par_beta;
	real<lower=0> s_alpha_hyper_par;
}

transformed data{
	vector[N_observed] n_m1_observed;
	vector[N_observed] inv_sqrt_n_observed;
	vector[N_observed] inv_scaled_s2_observed;

	n_m1_observed = n_observed - rep_vector(1, N_observed);
	inv_sqrt_n_observed = inv(sqrt(n_observed));
	inv_scaled_s2_observed = inv(s2_observed .* n_m1_observed);
}
			
parameters {
	real<lower=0> m_rho;
	real<lower=0> m_alpha;
	real<lower=0> s_rho;
	real<lower=0> s_alpha;
	real m_base;
	real s_base;
	vector[N_predict] m_f_tilde;
	vector[N_predict] s_f_tilde;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] m_predict;
	vector[N_predict] s_predict;
	// GP on mean of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha, m_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict = m_base + L_cov * m_f_tilde;
	// GP on log sigma of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha, s_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict = exp( s_base + L_cov * s_f_tilde);
}
			
model {
	m_rho ~ inv_gamma(m_rho_hyper_par_alpha, m_rho_hyper_par_beta);
	m_alpha ~ normal(0, m_alpha_hyper_par);	
	m_base ~ normal( 0 , 10 );
	s_rho ~ inv_gamma(s_rho_hyper_par_alpha, s_rho_hyper_par_beta);
	s_alpha ~ normal(0, s_alpha_hyper_par);	
	s_base ~ normal( 0 , 10 );	
	m_f_tilde ~ normal(0, 1);
	s_f_tilde ~ normal(0, 1);

	target+= normal_lpdf( y_observed | m_predict[observed_idx], s_predict[observed_idx] .* inv_sqrt_n_observed);
	target+= inv_chi_square_lpdf( s_predict[observed_idx] .* s_predict[observed_idx] .* inv_scaled_s2_observed | n_m1_observed );  	
}

generated quantities {
	vector[N_predict] exp_m_predict = exp( m_predict + s_predict .* s_predict .* rep_vector(0.5, N_predict) );  
}			

