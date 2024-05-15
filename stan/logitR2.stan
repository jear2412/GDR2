/* 
logit R2 

Prior over R2. Dirichlet decomposition of the total variance. 

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
  vector R2D2(vector z, vector sds_X, vector phi, real tau2) {
    /* Efficient computation of the R2D2 prior
    * Args:
    *   z: standardized population-level coefficients
    *   phi: local weight parameters
    *   tau2: global scale parameter (sigma is inside tau2)
    * Returns:
      *   population-level coefficients following the R2D2 prior
    */
    return z .* sqrt(phi * tau2) ./ sds_X;
  }
}
data {
  int<lower=1> N; // total number of observations
  vector[N] y; // response variable
  int<lower=1> p; // number of population-level effects (includes intercept!)
  matrix[N, p] X; // population-level design matrix (includes a column of ones)
  
  //real<lower=0> sigma;  // dispersion parameter
  
  // mean and covariance of logitphi
  vector[p - 2] mu_logitphi;
  cov_matrix[p - 2] sqcov_logiphi; //square root of the covariance matrix
  
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
  vector[N] yc;
  real ymean;
  
  matrix[Ntest, pc] Xctest; // centered version of X without an intercept
  vector[pc] means_Xtest; // column means of X before centering
  
  for (i in 2 : p) {
    means_X[i - 1] = mean(X[ : , i]);
    sds_X[i - 1] = sd(X[ : , i]);
    Xc[ : , i - 1] = X[ : , i] - means_X[i - 1];
  }
  
  ymean = mean(y);
  for (i in 1 : N) {
    yc[i] = y[i] - mean(y);
  }
  
  for (i in 2 : p) {
    means_Xtest[i - 1] = mean(Xtest[ : , i]);
    Xctest[ : , i - 1] = Xtest[ : , i] - means_Xtest[i - 1];
  }
}
parameters {
  // local parameters for the LNR2 prior
  vector[pc] zbeta;
  vector[pc - 1] zphi; //transform zphi to simplex 
  real Intercept; // temporary intercept for centered predictors
  // R2D2 shrinkage parameters
  real<lower=0, upper=1> R2LN_R2; // R2 parameter
  real<lower=0> sigma; // dispersion parameter
}
transformed parameters {
  vector[pc] beta; // population-level effects
  real R2_tau2; // global R2D2 scale parameter
  simplex[pc] R2_phi;
  vector[pc - 1] logitR2_phi;
  vector[pc] temp_logitR2D2phi;
  R2_tau2 = sigma ^ 2 * R2LN_R2 / (1 - R2LN_R2); //scaled by sigma
  // compute actual regression coefficients
  
  logitR2_phi = mu_logitphi + sqcov_logiphi * zphi; // There might be a more efficient way check, check brms 
  
  temp_logitR2D2phi[1] = 0;
  temp_logitR2D2phi[2 : pc] = logitR2_phi;
  R2_phi = softmax(temp_logitR2D2phi); // ?
  
  beta = R2D2(zbeta, sds_X, R2_phi, R2_tau2);
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma);
  }
  
  // priors including constants
  //R2 ~ Beta(R2mean, R2prec)
  target += beta_lpdf(R2LN_R2 | R2_mean * R2_prec, (1 - R2_mean)* R2_prec);
  
  target += std_normal_lpdf(zbeta); //normal distribution zbeta
  
  target += std_normal_lpdf(zphi); //normal distribution zphi
  
  target += student_t_lpdf(sigma | 3, 0, scale_sigma);
  
  
  target += normal_lpdf(Intercept | 0, 5); // Intercept  
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, beta);
  
  vector[N] log_lik;
  array[N] real y_tilde;
  vector[N] mu_tilde = rep_vector(0.0, N) + ymean + Intercept + Xc * beta;
  
  vector[Ntest] log_lik_test;
  array[Ntest] real y_tilde_test;
  vector[Ntest] mu_tilde_test = rep_vector(0.0, Ntest) + b_Intercept + Xtest[,2:p] * beta;
  
  // R2
  real<lower=0, upper=1> pred_R2 = variance(mu_tilde) / (variance(mu_tilde) + sigma^2 );
  real<lower=0, upper=1> pred_R2_test = variance(mu_tilde_test) / (variance(mu_tilde_test) + sigma^2 );
  
  real tau=1.0;
  
  // lambdas
  
  vector[pc] lambda = R2_phi * R2_tau2; //variances of betas
  
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


