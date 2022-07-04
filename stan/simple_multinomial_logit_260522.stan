// Example at https://mc-stan.org/docs/2_29/stan-users-guide/multi-logit.html
// https://eleafeit.com/posts/2021-05-23-parameterization-of-multinomial-logit-models-in-stan/

functions
{
}

data
{
  int<lower=1> N; // # of participants
  int<lower=2> K; // # of classes
  int<lower=1> D; // # of predictors
  array[N] int<lower=1, upper=K> y;     // label for each participant
  matrix[N, D] x; // predictors
}

transformed data
{
}

parameters
{
  matrix[D, K-1] beta_raw;
}

transformed parameters
{
  vector[D] zeros = rep_vector(0, D);
  matrix[D, K] beta = append_col(beta_raw, zeros);
}

model
{
  matrix[N, K] x_beta= x*beta;
  to_vector(beta_raw) ~ normal(0,5);

  for (n in 1:N)
  {
    y[n] ~ categorical_logit(x_beta[n]'); // ' is used as Traspose
  }
}

generated quantities
{
  // Probabilities of belonging in different classes
  // As a function of LALALA
  softmax((1,0) * beta)
  softmax((1,1) * beta)
}
