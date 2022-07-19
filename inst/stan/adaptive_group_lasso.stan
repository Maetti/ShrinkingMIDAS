//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
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
  vector[2] pr_lambda;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[nZ] theta;
  real<lower=0> sigma_2;
  real<lower=0> tau[nG];
  real<lower=0> lambda[nG];
}

transformed parameters {

  vector[nZ] tr_Sigma;
  real<lower=0> lambda_half[nG];

  for (i in 1:nZ) {
    tr_Sigma[i] = sigma_2 * square(tau[gInd[i]]);
  };

  for (i in 1:nG) {
    lambda_half[i] = lambda[i] / 2;
  };

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  for (i in 1:nG) {
    lambda[i] ~ gamma(pr_lambda[1], pr_lambda[2]);
  };

  // for (i in 1:nG) {
  //  tau[i] ~ gamma((gSize[i] + 1)/2, lambda[i]/2);
  // };

  tau ~ gamma((gSize + 1)/2, lambda_half);

  sigma_2 ~ inv_gamma(pr_sigma[1], pr_sigma[2]);

  theta ~ normal(0, tr_Sigma);

  y_train ~ normal(x_train * theta, sigma_2);

}


generated quantities {
  real y_rep[nY_train];
  real y_pred[nY_test];
  vector[nY_test] log_lik;
  for (i in 1:nY_test) {
    y_pred[i] = normal_rng(x_test[i] * theta, sigma_2);
    log_lik[i] = normal_lpdf(y_test[i] | x_test[i] * theta, sigma_2);
  };

  for (j in 1:nY_train) {
    y_rep[j] = normal_rng(x_train[j] * theta, sigma_2);
  };
}
