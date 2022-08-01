
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
	real baseline;
	real sex;
	real loc;
	vector[AGE_N] phi; 		
	real<lower=0> sigma;
}
			
transformed parameters{
	vector[N] p_logit;
	for ( i in 1:N ) 
	{
		p_logit[i] = baseline + sex * SEX[i] + loc * LOC[i] + sigma * phi[AGE[i]];
	}
}
			
model{	
	baseline ~ normal( 0 , 10 );
	loc ~ normal( 0 , 3 );
	sex ~ normal( 0 , 3 );	
	phi ~ icar_normal_lpdf(AGE_N , node1 , node2);
	sigma ~ cauchy(0, 1);
	K ~ binomial_logit( TOTAL , p_logit );
}
			
generated quantities{
	vector[N] p;
	p= inv_logit(p_logit);
}
