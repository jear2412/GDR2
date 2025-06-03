/* /Normal Gamma prior Griffin,Brown 2010

see section 3.3 in the paper
beta_i | psi2_i ~ N(0, psi2_i ) (psi2 is a variance here)
psi2_i | lambda, gamma ~   Gamma(c, d)
c= lambda
d= 1/(2*gamma^2)

nu = 2 lambda gamma^2  ~ Inverse Gamma(2, M)
lambda ~ exp(1)

gamma | lambda = sqrt(nu/(2.0*lambda))

Cases: 
1) X nonsingular 
M = || ols(beta) ||_2^2 
2) X singular
M is = 1/n || MMLS(beta) ||_2^2 , MMLS : Minimum Length Least squares

M can also be an educated guess about the l2 norm of beta. 

We diverge in the following:
1. Sigma They propose an improper prior 
2.Intercept Improper prior

*/
data {
  int<lower=1> N; // Number of observations
  int<lower=1> p; // Number of covariates (includes intercept)
  matrix[N, p] X; // Includes a column of 1s for intercept
  vector[N] y;
  //real<lower=0> sigma; // value for sigma
  real<lower=0> M;// value for M: educated guess about l2 norm of beta
  
  
  //---- test data
  int<lower=1> Ntest; // total number of observations
  vector[Ntest] ytest; // test set
  matrix[Ntest, p] Xtest; // population-level design matrix including column of 1s
  real<lower=0> scale_sigma; // scale for the distribution of sigma
  
  int prior_only;  // should the likelihood be ignored?
  
}

transformed data {
  int pc = p - 1;
  matrix[N, pc] Xc; // centered version of X without an intercept
  vector[pc] means_X; // column means of X before centering
  vector[pc] sds_X; // sds of X before centering
  vector[pc] inv_sds_X;
  vector[N] yc;
  real ymean;
  
  for (i in 2 : p) {
    means_X[i - 1] = mean(X[ : , i]);
    sds_X[i - 1] = sd(X[ : , i]);
    inv_sds_X[i-1] = 1.0 / sds_X[i-1];
    Xc[ : , i - 1] = X[ : , i] - means_X[i - 1];
  }
  
  ymean = mean(y);
  for (i in 1 : N) {
    yc[i] = y[i] - ymean;
  }
  
 
}

parameters {
  real Intercept;  // temporary intercept for centered predictors
  vector[pc] zbeta;
  vector<lower=0>[pc] psi2; // variances of b
  real<lower=0> lambda;
  real<lower=0> nu;
  real<lower=0> sigma;
}

transformed parameters {
  vector[pc] beta;
  real<lower=0> gamma;
  
  gamma= sqrt(nu/(2.0*lambda)); // nu = 2 *lambda*gamma^2
  
  if(!prior_only) {
    beta= sigma* zbeta .* sqrt(psi2) .* inv_sds_X; 
  }else{
    beta= sigma* zbeta .* sqrt(psi2)  ;
  }
  
}

model {
  
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma); //Intercept+Xc*beta

  }
  
  // priors including constants
  target += std_normal_lpdf(zbeta); // prior for zb
 
  target += exponential_lpdf( lambda | 1);  // prior for lambda Exponential(1) 
  
  target+= inv_gamma_lpdf(nu | 2, M); // prior for nu Inverse Gamma(2, M), expected value is M , relationship bw nu and gamma 

  target += gamma_lpdf( psi2 | lambda, 1/(2.0*pow(gamma,2))  ); // psi2 ~ Gamma(lambda, 1 / (2 * gamma^2))  => E[psi2] = lambda * 2 * gamma^2
  
  target += student_t_lpdf(sigma | 3, 0, scale_sigma );
    
  
  target += normal_lpdf(Intercept | 0, 5);  // Intercept
  
}

generated quantities {
  // actual population-level intercept
  real b_Intercept = ymean + Intercept - dot_product(means_X, beta);

  vector[N] log_lik;
  array[N] real y_tilde;
  vector[N] mu_tilde = rep_vector(0.0, N) +  b_Intercept + X[, 2:p] * beta;

  vector[Ntest] log_lik_test;
  array[Ntest] real y_tilde_test;
  vector[Ntest] mu_tilde_test = rep_vector(0.0, Ntest) + b_Intercept + Xtest[,2:p] * beta;
  

  //Predictive R2
  real<lower=0, upper=1> pred_R2 = variance(mu_tilde) / (variance(mu_tilde) + sigma^2 );
  real<lower=0, upper=1> pred_R2_test = variance(mu_tilde_test) / (variance(mu_tilde_test) + sigma^2 );
 
  //---y_tilde calc
  for (n in 1 : N) {
    log_lik[n] = normal_lpdf(y[n] | mu_tilde[n], sigma);
    y_tilde[n] = normal_rng(mu_tilde[n], sigma); //copy and paste model (executed once per sample) 
  }
  
  //---y_tilde test calc
  for (n in 1 : Ntest) {
    log_lik_test[n] = normal_lpdf(ytest[n] | mu_tilde_test[n], sigma);
    y_tilde_test[n] = normal_rng(mu_tilde_test[n], sigma); //copy and paste model (executed once per sample) 
  }
  
  
}

