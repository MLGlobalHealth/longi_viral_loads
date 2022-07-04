// Explain features of the model here:
functions
{
  // Annoying that I cannot add int to integer array but loop saves me here
  array[] int add2array( array[] int arr, int i, int L ) {
    array[L] int out;
    for (idx in 1:L) {
      out[idx] = arr[idx] + i;
    }
    return out;
  }
  
}

data
{
  // int<lower=1> NF; // # of Female Participants with at least one visit pair
  // int<lower=1> NM; // # of Female Participants with at least one visit pair
  int<lower=1> NP; // # of visitpairs
  array[NP] int<lower=0, upper=1> v_final;
  array[NP] int<lower=0, upper=1> v_start;
  array[NP] int<lower=0, upper=1> zF;
  // alternatively, report the ID for the subject, and then use look-up tables
  // to get different info when needed.
}

transformed data
{
  array[NP] int idx_sex;
  array[NP] int idx_start;

  idx_sex = add2array(zF, 1, NP); 
  idx_start = add2array(v_start, 1, NP); 
}

parameters
{
  matrix<lower=0, upper=1>[2,2]  pi; // Nrow = Number of groups. Ncol = possible v_starttart
}

transformed parameters
{
  vector[NP] mu;  // Means to be passed to the bernoulli pdf

  for (idx in 1:NP) {
    mu[idx] = pi[idx_sex[idx], idx_start[idx] ];
  }

}

model
{
  v_final ~ bernoulli(mu);
}

generated quantities
{
} 
