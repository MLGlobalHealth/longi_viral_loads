
functions{			
	real icar_normal_lpdf( vector phi, int N, int[] node1, int[] node2)
	{			
		return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf( sum(phi) | 0, 0.001 * N);
	}
}
		
data{
	int<lower=1> N; 
	int<lower=1> TOTAL[N];    
	int<lower=0> K[N];
	int<lower=0> AGE_N;  
	int<lower=0, upper=AGE_N> AGE[N];  
	int<lower=0,upper=1> SEX[N];
	int<lower=0,upper=1> LOC[N];	 
	int<lower=0> N_edges; 						// number of related age groups
	int<lower=1, upper=AGE_N> node1[N_edges];  	// node1[i] adjacent to node2[i]
	int<lower=1, upper=AGE_N> node2[N_edges];  	// and node1[i] < node2[i]
}
		
parameters{
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;
	vector[AGE_N] phi_sex0_loc0;
	vector[AGE_N] phi_sex1_loc0;
	vector[AGE_N] phi_sex0_loc1;
	vector[AGE_N] phi_sex1_loc1; 			 		
	real<lower=0> sigma_loc0;
	real<lower=0> sigma_loc1;
}
		
transformed parameters{
	vector[N] p_logit;
	for ( i in 1:N ) 
	{
		p_logit[i] = sex0_loc0 * (1-SEX[i]) * (1-LOC[i]) + 
					 sex1_loc0 * SEX[i] * (1-LOC[i]) + 
					 sex0_loc1 * (1-SEX[i]) * LOC[i] + 
					 sex1_loc1 * SEX[i] * LOC[i] +
						sigma_loc0 * phi_sex0_loc0[AGE[i]] * (1-SEX[i]) * (1-LOC[i]) +
						sigma_loc0 * phi_sex1_loc0[AGE[i]] * SEX[i] * (1-LOC[i]) +
						sigma_loc1 * phi_sex0_loc1[AGE[i]] * (1-SEX[i]) * LOC[i] +
						sigma_loc1 * phi_sex1_loc1[AGE[i]] * SEX[i] * LOC[i];
	}
}
		
model{	
	sex0_loc0 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
	phi_sex0_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex0_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	sigma_loc0 ~ normal(0.0, 1);
	sigma_loc1 ~ normal(0.0, 1);
	//sigma_loc0 ~ cauchy(0, 1);
	//sigma_loc1 ~ cauchy(0, 1);
	K ~ binomial_logit( TOTAL , p_logit );
}
		
generated quantities{
	vector[N] p;
	p= inv_logit(p_logit);
}

