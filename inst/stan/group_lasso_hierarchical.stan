// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// data
data {
  int<lower=0> nY_train;
  int<lower=0> nY_test;
  int<lower=0> nZ;
  int<lower=0> nG;
  real y_train[nY_train];
  real y_test[nY_test];
  matrix[nY_train, nZ] x_train;
  matrix[nY_test, nZ] x_test;
  vector[nG] gSize;
  int gInd[nZ];
  vector[2] pr_sigma;
}


// The parameters accepted by the model
parameters {
  vector[nZ] theta;
  real<lower=0> sigma2;

  // shrinkage parameters
  real<lower=0> tau; // global
  real<lower=0> delta[nG]; // group
  real<lower=0> lambda[nZ]; // local
}

transformed parameters {
  real<lower=0> sigma;
  vector[nZ] tr_Sigma;

  sigma = sqrt(sigma2);

  for (i in 1:nZ) {
    tr_Sigma[i] = sigma2 * square(lambda[i]) * square(delta[gInd[i]]);
  };

}

// The model to be estimated.
model {

  // global
  tau ~ cauchy(0, 1);

  // group
  delta ~ exponential(1);

  // local
  lambda ~ exponential(1);

  // sigma
  sigma2 ~ inv_gamma(pr_sigma[1], pr_sigma[2]);

  // coefficient
  theta ~ normal(0, tr_Sigma);

  // model
  y_train ~ normal(x_train * theta, sigma);

}


generated quantities {
  // posterior predictive checking (how good does the model fit the train data)
  real y_rep[nY_train] = normal_rng(x_train * theta, sigma);
  
  // posterior hold-out predictive sampling (fit to new data) and log density
  real y_pred[nY_test] = normal_rng(x_test * theta, sigma);
  real log_lik = normal_lpdf(y_test | x_test * theta, sigma);
}
