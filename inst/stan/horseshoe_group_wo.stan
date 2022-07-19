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

// The parameters accepted by the model.
parameters {
  vector[nZ] theta;
  real<lower=0> sigma_2;

  // shrinkage parameters
  real<lower=0> tau; // global
  real<lower=0> delta[nG]; // group
  real<lower=0> lambda[nZ]; // local
}

transformed parameters {
  real<lower=0> sigma;
  vector[nZ] tr_Sigma;

  sigma = sqrt(sigma_2);

  for (i in 1:nZ) {
    tr_Sigma[i] = sigma_2 * lambda[i] * delta[gInd[i]];
  };

}

// The model to be estimated.
model {

  // global
  tau ~ cauchy(0, 1);

  // group
  for (i in 1:nG) {
    delta[i] ~ cauchy(0, 1);
  };

  // local
  for (i in 1:nZ) {
    lambda[i] ~ cauchy(0, 1);
  };


  sigma_2 ~ inv_gamma(pr_sigma[1], pr_sigma[2]);

  theta ~ normal(0, tr_Sigma);

  y_train ~ normal(x_train * theta, sigma);

}

