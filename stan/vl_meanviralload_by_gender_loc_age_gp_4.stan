
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed_00;
	int<lower=1> N_observed_10;
	int<lower=1> N_observed_01;
	int<lower=1> N_observed_11;
	int<lower=1, upper=N_predict> observed_idx_00[N_observed_00];
	int<lower=1, upper=N_predict> observed_idx_10[N_observed_10];
	int<lower=1, upper=N_predict> observed_idx_01[N_observed_01];
	int<lower=1, upper=N_predict> observed_idx_11[N_observed_11];
	vector<lower=0>[N_observed_00] y_observed_00;
	vector<lower=0>[N_observed_10] y_observed_10;
	vector<lower=0>[N_observed_01] y_observed_01;
	vector<lower=0>[N_observed_11] y_observed_11;
	vector<lower=0>[N_observed_00] s2_observed_00; // with denominator (n-1)
	vector<lower=0>[N_observed_10] s2_observed_10; // with denominator (n-1)
	vector<lower=0>[N_observed_01] s2_observed_01; // with denominator (n-1)
	vector<lower=0>[N_observed_11] s2_observed_11; // with denominator (n-1)
	vector<lower=1>[N_observed_00] n_observed_00;
	vector<lower=1>[N_observed_10] n_observed_10;
	vector<lower=1>[N_observed_01] n_observed_01;
	vector<lower=1>[N_observed_11] n_observed_11;
	real<lower=0> m_rho_hyper_par_00;
	real<lower=0> m_alpha_hyper_par_00;
	real<lower=0> s_rho_hyper_par_00;
	real<lower=0> s_alpha_hyper_par_00;
	real<lower=0> m_rho_hyper_par_10;
	real<lower=0> m_alpha_hyper_par_10;
	real<lower=0> s_rho_hyper_par_10;
	real<lower=0> s_alpha_hyper_par_10;
	real<lower=0> m_rho_hyper_par_01;
	real<lower=0> m_alpha_hyper_par_01;
	real<lower=0> s_rho_hyper_par_01;
	real<lower=0> s_alpha_hyper_par_01;
	real<lower=0> m_rho_hyper_par_11;
	real<lower=0> m_alpha_hyper_par_11;
	real<lower=0> s_rho_hyper_par_11;
	real<lower=0> s_alpha_hyper_par_11;
}

transformed data{
	vector[N_observed_00] n_m1_observed_00;
	vector[N_observed_00] inv_sqrt_n_observed_00;
	vector[N_observed_00] inv_scaled_s2_observed_00;
	vector[N_observed_10] n_m1_observed_10;
	vector[N_observed_10] inv_sqrt_n_observed_10;
	vector[N_observed_10] inv_scaled_s2_observed_10;
	vector[N_observed_01] n_m1_observed_01;
	vector[N_observed_01] inv_sqrt_n_observed_01;
	vector[N_observed_01] inv_scaled_s2_observed_01;
	vector[N_observed_11] n_m1_observed_11;
	vector[N_observed_11] inv_sqrt_n_observed_11;
	vector[N_observed_11] inv_scaled_s2_observed_11;

	n_m1_observed_00 = n_observed_00 - rep_vector(1, N_observed_00);
	inv_sqrt_n_observed_00 = inv(sqrt(n_observed_00));
	inv_scaled_s2_observed_00 = inv(s2_observed_00 .* n_m1_observed_00);
	n_m1_observed_10 = n_observed_10 - rep_vector(1, N_observed_10);
	inv_sqrt_n_observed_10 = inv(sqrt(n_observed_10));
	inv_scaled_s2_observed_10 = inv(s2_observed_10 .* n_m1_observed_10);
	n_m1_observed_01 = n_observed_01 - rep_vector(1, N_observed_01);
	inv_sqrt_n_observed_01 = inv(sqrt(n_observed_01));
	inv_scaled_s2_observed_01 = inv(s2_observed_01 .* n_m1_observed_01);
	n_m1_observed_11 = n_observed_11 - rep_vector(1, N_observed_11);
	inv_sqrt_n_observed_11 = inv(sqrt(n_observed_11));
	inv_scaled_s2_observed_11 = inv(s2_observed_11 .* n_m1_observed_11);
}
			
parameters {
	real<lower=0> m_rho_00;
	real<lower=0> m_alpha_00;
	real<lower=0> s_rho_00;
	real<lower=0> s_alpha_00;
	real m_base_00;
	real s_base_00;
	vector[N_predict] m_f_tilde_00;
	vector[N_predict] s_f_tilde_00;	
	real<lower=0> m_rho_10;
	real<lower=0> m_alpha_10;
	real<lower=0> s_rho_10;
	real<lower=0> s_alpha_10;
	real m_base_10;
	real s_base_10;
	vector[N_predict] m_f_tilde_10;
	vector[N_predict] s_f_tilde_10;	
	real<lower=0> m_rho_01;
	real<lower=0> m_alpha_01;
	real<lower=0> s_rho_01;
	real<lower=0> s_alpha_01;
	real m_base_01;
	real s_base_01;
	vector[N_predict] m_f_tilde_01;
	vector[N_predict] s_f_tilde_01;	
	real<lower=0> m_rho_11;
	real<lower=0> m_alpha_11;
	real<lower=0> s_rho_11;
	real<lower=0> s_alpha_11;
	real m_base_11;
	real s_base_11;
	vector[N_predict] m_f_tilde_11;
	vector[N_predict] s_f_tilde_11;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] m_predict_00;
	vector[N_predict] s_predict_00;
	vector[N_predict] m_predict_10;
	vector[N_predict] s_predict_10;
	vector[N_predict] m_predict_01;
	vector[N_predict] s_predict_01;
	vector[N_predict] m_predict_11;
	vector[N_predict] s_predict_11;
	// GP on mean of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_00, m_rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_00 = m_base_00 + L_cov * m_f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_10, m_rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_10 = m_base_10 + L_cov * m_f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_01, m_rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_01 = m_base_01 + L_cov * m_f_tilde_01;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_11, m_rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_11 = m_base_11 + L_cov * m_f_tilde_11;
	// GP on log sigma of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_00, s_rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_00 = exp( s_base_00 + L_cov * s_f_tilde_00);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_10, s_rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_10 = exp( s_base_10 + L_cov * s_f_tilde_10);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_01, s_rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_01 = exp( s_base_01 + L_cov * s_f_tilde_01);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_11, s_rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_11 = exp( s_base_11 + L_cov * s_f_tilde_11);
}
			
model {
	m_rho_00 ~ normal(0, m_rho_hyper_par_00);
	m_alpha_00 ~ normal(0, m_alpha_hyper_par_00);	
	m_base_00 ~ normal( 0 , 10 );
	s_rho_00 ~ normal(0, s_rho_hyper_par_00);
	s_alpha_00 ~ normal(0, s_alpha_hyper_par_00);	
	s_base_00 ~ normal( 0 , 10 );	
	m_f_tilde_00 ~ normal(0, 1);
	s_f_tilde_00 ~ normal(0, 1);

	m_rho_10 ~ normal(0, m_rho_hyper_par_10);
	m_alpha_10 ~ normal(0, m_alpha_hyper_par_10);	
	m_base_10 ~ normal( 0 , 10 );
	s_rho_10 ~ normal(0, s_rho_hyper_par_10);
	s_alpha_10 ~ normal(0, s_alpha_hyper_par_10);	
	s_base_10 ~ normal( 0 , 10 );	
	m_f_tilde_10 ~ normal(0, 1);
	s_f_tilde_10 ~ normal(0, 1);

	m_rho_01 ~ normal(0, m_rho_hyper_par_01);
	m_alpha_01 ~ normal(0, m_alpha_hyper_par_01);	
	m_base_01 ~ normal( 0 , 10 );
	s_rho_01 ~ normal(0, s_rho_hyper_par_01);
	s_alpha_01 ~ normal(0, s_alpha_hyper_par_01);	
	s_base_01 ~ normal( 0 , 10 );	
	m_f_tilde_01 ~ normal(0, 1);
	s_f_tilde_01 ~ normal(0, 1);

	m_rho_11 ~ normal(0, m_rho_hyper_par_11);
	m_alpha_11 ~ normal(0, m_alpha_hyper_par_11);	
	m_base_11 ~ normal( 0 , 10 );
	s_rho_11 ~ normal(0, s_rho_hyper_par_11);
	s_alpha_11 ~ normal(0, s_alpha_hyper_par_11);	
	s_base_11 ~ normal( 0 , 10 );	
	m_f_tilde_11 ~ normal(0, 1);
	s_f_tilde_11 ~ normal(0, 1);

	target+= normal_lpdf( y_observed_00 | m_predict_00[observed_idx_00], s_predict_00[observed_idx_00] .* inv_sqrt_n_observed_00);
	target+= inv_chi_square_lpdf( s_predict_00[observed_idx_00] .* s_predict_00[observed_idx_00] .* inv_scaled_s2_observed_00 | n_m1_observed_00 );
	target+= normal_lpdf( y_observed_10 | m_predict_10[observed_idx_10], s_predict_10[observed_idx_10] .* inv_sqrt_n_observed_10);
	target+= inv_chi_square_lpdf( s_predict_10[observed_idx_10] .* s_predict_10[observed_idx_10] .* inv_scaled_s2_observed_10 | n_m1_observed_10 );
	target+= normal_lpdf( y_observed_01 | m_predict_01[observed_idx_01], s_predict_01[observed_idx_01] .* inv_sqrt_n_observed_01);
	target+= inv_chi_square_lpdf( s_predict_01[observed_idx_01] .* s_predict_01[observed_idx_01] .* inv_scaled_s2_observed_01 | n_m1_observed_01 );
	target+= normal_lpdf( y_observed_11 | m_predict_11[observed_idx_11], s_predict_11[observed_idx_11] .* inv_sqrt_n_observed_11);
	target+= inv_chi_square_lpdf( s_predict_11[observed_idx_11] .* s_predict_11[observed_idx_11] .* inv_scaled_s2_observed_11 | n_m1_observed_11 );  	
}

generated quantities {
	vector[N_predict] exp_m_predict_00;
	vector[N_predict] exp_m_predict_10;
	vector[N_predict] exp_m_predict_01;
	vector[N_predict] exp_m_predict_11; 
	exp_m_predict_00 = exp( m_predict_00 + s_predict_00 .* s_predict_00 .* rep_vector(0.5, N_predict) );
	exp_m_predict_10 = exp( m_predict_10 + s_predict_10 .* s_predict_10 .* rep_vector(0.5, N_predict) );
	exp_m_predict_01 = exp( m_predict_01 + s_predict_01 .* s_predict_01 .* rep_vector(0.5, N_predict) );
	exp_m_predict_11 = exp( m_predict_11 + s_predict_11 .* s_predict_11 .* rep_vector(0.5, N_predict) );  
}	
