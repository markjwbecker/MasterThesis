rm(list = ls())
set.seed(123)
setwd("C:/Users/markj/Documents/R files")
library(readxl)
library(rstan)
library(MTS)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

Phi1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
Sigma=matrix(c(4,0.8,0.8,1),2,2)
Theta1=matrix(c(-0.5,0,0,-0.6),2,2)
sim_data=VARMAsim(250,arlags=c(1),malags=c(1),phi=Phi1,theta=Theta1,sigma=Sigma)

z=sim_data$series
plot.ts(z)

stan_code <- "

data {
  int<lower=1> T; //number of observations
  int<lower=2> m; //number of variables
  matrix[T, m] z; //data
}


parameters {
  matrix[m, m] Phi;
  matrix[m, m] Theta;
  cov_matrix[m] Sigma;
}

model {

  Sigma ~ inv_wishart(m+2, diag_matrix(rep_vector(1, m)));
  to_vector(Phi) ~ multi_normal(rep_vector(0, m*m), diag_matrix(rep_vector(1, m*m)));
  to_vector(Theta) ~ multi_normal(rep_vector(0, m*m), diag_matrix(rep_vector(1, m*m)));
  
  matrix[T, m] nu;         // prediction for time t
  matrix[T, m] err;        // error for time t
  
  nu[1] = rep_row_vector(0, m);      // assume err[0] == 0, mu=0
  err[1] = z[1] - nu[1];
  
  for (t in 2:T) {
    nu[t] = (Phi * z[t - 1]' - Theta * err[t - 1]')';  // Fix type mismatch
    err[t] = z[t] - nu[t];
}

  for (t in 1:T) {
    err[t] ~ multi_normal(rep_vector(0, m), Sigma);
}
  
}
"

data <- list(z=z,m=ncol(z),T=nrow(z))

n_chains <- 2
iter <- 5000
warmup <- 1500
fit <- stan(model_code=stan_code,data=data,chains=n_chains,iter=iter,warmup=warmup, verbose = TRUE)




Phi1_est <- matrix(get_posterior_mean(fit)[1:4,3],2,2,byrow=TRUE)
Theta1_est <- matrix(get_posterior_mean(fit)[5:8,3],2,2,byrow=TRUE)
Sigma1_est <- matrix(get_posterior_mean(fit)[9:12,3],2,2,byrow=TRUE)


Phi1
Phi1_est

Theta1
Theta1_est

Sigma
Sigma_est



