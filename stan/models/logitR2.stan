/* 
logit R2 

Prior over R2. logit decomposition of the total variance. 
Proportion of explained variance R2
R2 ~ Beta(R2mean, R2prec)
The usual parametrization of the beta distribution is recovered by setting
a1= R2mean*R2prec
a2= (1-R2mean)*R2prec

Explained variance tau2
tau2 = R2/(1-R2) ~ BetaPrime(R2mean, R2prec)
phi ~ logit normal (mu, sigma)
lambda^2 = phi* tau^2 Proportion of explained variance

beta ~ Normal(0, sigma^2*phi*tau^2 )

Prior for sigma
Half Student t is recommended 

*/
functions {
  /* Efficient computation of the R2 prior
   * Args:
   *   z: standardized population-level coefficients
   *   phi: local weight parameters
   *   tau2: global scale parameter
   * Returns:
   *   population-level coefficients following the R2 prior
   */
  vector R2D(vector z, vector inv_sds_X, vector phi, real tau2, int prior_only) {
    if(!prior_only){
      return z .* sqrt(phi * tau2) .* inv_sds_X;   // scale by the sds of X
    } else {
      return z .* sqrt(phi * tau2);  // no scales needed
    }
    
  }
}
data {
  int<lower=1> N; // total number of observations
  vector[N] y; // response variable
  int<lower=1> p; // number of population-level effects (includes intercept!)
  matrix[N, p] X; // population-level design matrix (includes a column of ones)
  vector[p -1] mu_logitphi; // mean and covariance of logitphi
  cholesky_factor_cov[p - 1] sigma_logitphi; //square root of the covariance matrix
  
  //---- test data
  int<lower=1> Ntest; // total number of observations
  vector[Ntest] ytest; // test set
  matrix[Ntest, p] Xtest; // population-level design matrix including column of 1s
  
  // data for the R2 prior
  real<lower=0> R2_mean; // mean of the R2 prior
  real<lower=0> R2_prec; // precision of the R2 prior
  
  real<lower=0> scale_sigma; // scale for the distribution of sigma
  int prior_only; // should the likelihood be ignored?
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
  // local parameters for the LNR2 prior
  vector[pc] zbeta;
  vector[pc] zphi; //transform zphi to simplex 
  real Intercept; // temporary intercept for centered predictors
  
  // GDR2 shrinkage parameters
  real<lower=0, upper=1> R2LN_R2; // R2 parameter
  real<lower=0> sigma; // dispersion parameter
}

transformed parameters {
  vector[pc] beta; // population-level effects
  real R2_tau2; // global R2D2 scale parameter
  simplex[pc] R2_phi;
  vector[pc] logitR2_phi;
  
  R2_tau2 = sigma ^ 2 * R2LN_R2 / (1 - R2LN_R2); //scaled by sigma
  // compute actual regression coefficients
  
  // transform zphi 
  logitR2_phi = mu_logitphi + sigma_logitphi * zphi; 
  
  R2_phi = softmax(logitR2_phi); 

  beta = R2D(zbeta, inv_sds_X, R2_phi, R2_tau2, prior_only); //apply R2 decomposition
}
model {
  
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma);
  }
  

  target += beta_lpdf(R2LN_R2 | R2_mean * R2_prec, (1 - R2_mean)* R2_prec); //beta distribution for R2LN_R2
  
  target += std_normal_lpdf(zbeta); //normal distribution zbeta
  
  target += std_normal_lpdf(zphi); //normal distribution zphi
  
  target += student_t_lpdf(sigma | 3, 0, scale_sigma);
  
  target += normal_lpdf(Intercept | 0, 5); // Intercept  
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
  
  // R2
  real<lower=0, upper=1> pred_R2 = variance(mu_tilde) / (variance(mu_tilde) + sigma^2 );
  real<lower=0, upper=1> pred_R2_test = variance(mu_tilde_test) / (variance(mu_tilde_test) + sigma^2 );
  
  // lambdas
  vector[pc] lambda = R2_phi * R2_tau2; //variances of betas
  
  //--- y_tilde calc
  for (n in 1 : N) {
    log_lik[n] = normal_lpdf(y[n] | mu_tilde[n], sigma);
    y_tilde[n] = normal_rng(mu_tilde[n], sigma); 
  }
  
  //---y_tilde test calc
  for (n in 1 : Ntest) {
    log_lik_test[n] = normal_lpdf(ytest[n] | mu_tilde_test[n], sigma);
    y_tilde_test[n] = normal_rng(mu_tilde_test[n], sigma); 
  }
  
}


