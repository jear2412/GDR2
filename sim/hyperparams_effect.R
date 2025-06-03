library(brms)
library(bayesplot)
library(cmdstanr)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(tidyverse)
library(ggcorrplot)
library(foreach)
library(future)
library(ggridges)
library(Matrix)
library(doParallel)
library(doFuture)
library(doRNG)

source("scripts/aux_functions/all_auxfunctions.R") #script containing auxiliary functions
source("scripts/aux_functions/coef_functions.R") #script containing coefficients auxiliary functions
source("scripts/aux_functions/cov_functions.R") #script containing covariance auxiliary functions
source("scripts/aux_functions/dgp_functions.R") #script containing dgp auxiliary functions
source("scripts/aux_functions/summary_functions.R") #script containing summary auxiliary functions

# Design -----

# normal means problem

n= 50 #size of training data
ntest= n # size of test data 
ps= 200 #num overall of coefficients 
rhos=c(0) #correlation of X, X ~ MVN
nus= c(0.1) #sparsity level 
type=c("ar") #type of correlation matrix used
alpha= 0 #real intercept
R2= c(0.80)  # signal to noise ratio
dgp = c("normal_means_dgp")
omega_flag = c("true")

# functions to generate coefficients
gencoef_function = as.character(c("gen_coef_normal_means" ))

#simulation conditions
sim_conds <- expand.grid(n=n, ntest=ntest, p=ps, 
                         rho=rhos, nu=nus, 
                         type= as.character(type), 
                         alpha= alpha, R2= R2, 
                         dgp= dgp, 
                         gencoef_function = gencoef_function,
                         omega_flag = omega_flag,
                         proj_pred_flag = 0, 
                         seed = 343)

sim_cond <- sim_conds[1,]
#params <- normal_means_dgp( params = list(sim_cond = sim_cond) )
params <- general_dgp( params = list(sim_cond = sim_cond) )

niters = 2000
prior_flag = 1

# Create empty list to store results
results <- list()

sigma_ln <- sigma_ln <- c(1e-3, 1e-2, 1e-1, 1, 10, 100)

#------ Model 

file <- file.path("stan/models", "logitR2.stan")

mod <- cmdstan_model(file)

# Use foreach to run the models in parallel
# Set up parallel backend
registerDoFuture()
plan(multisession, workers = 12)  # Adjust the number of workers as needed

results <- foreach(sigma_val = sigma_ln, .packages = c("cmdstanr", "dplyr", "brms")) %dopar% {
  p = params$p
  d = p #dimension
  # Modify the sigma_logitphi for each sigma_ln value
  
  scale= sqrt(log(p)/p)
  scale = 1
  
  sigma_val = sigma_val * scale 
  
  rho <- 0
  blockA <- diag(sigma_val, 5 ) %*%( (1-rho)*diag(1, 5)+rho*matrix(1, 5, 5)) %*% diag(sigma_val, 5 ) 
  blockB <- diag(sigma_val^2, d-5)
  sigma_logitphi <- bdiag(blockA, blockB)
  
  #sigma_logitphi <- diag(1, p)
  #sigma_logitphi <- diag(sigma_val, p ) %*% sigma_logitphi %*% diag(sigma_val, p)
  chol_sigma_logitphi <- as.matrix(t(chol(sigma_logitphi)))
  
  mu_logit = c(rep(0, d)) 
  #mu_logit = c(c(2,2,2), c(1,1), rep(0, d-5)) 
  
  # Prepare the data list for Stan
  dat <- list(
    N= params$n,
    p= params$p+1,
    X=cbind(rep(1, n), params$X),
    y=as.numeric(params$y),
    Ntest=params$ntest, 
    Xtest=cbind(rep(1, params$ntest), params$Xtest),
    ytest=as.numeric(params$ytest),
    scale_sigma= sd(params$y),
    mu_logitphi= mu_logit,
    chol_sigma_logitphi= chol_sigma_logitphi,
    R2_mean= 0.5, 
    R2_prec=1,
    R2_alpha = rep(0.5, p),
    eta=1, 
    prior_only=prior_flag
  )
  
  # Fit the model using cmdstanr
  
  fit <- mod$sample(
    data = dat,
    chains = 1,
    refresh = 250,
    seed = sim_cond$seed,
    #parallel_chains = 4,
    adapt_delta = 0.95,
    iter_sampling = niters
  )
  
  sm <- fit$draws(format = "df")
  
  sm <- sm %>% 
        mutate( sigma_ln = sigma_val)
  
  # Return the results
  sm
}

# Points to discuss:
# Scenarios: changing mu and sigma. I can include a block of 5 x 5 correlated.
# When mu = 0  the posterior of beta1 is shrunken heavily as sigma_ln increases
# I give another structur to mu as well (discuss) 

# Case 0: prior draws 
# Case 1: n >> p (200 vs 50)

# Case 2: n << p (50 vs 200) 
# As sigma is lower the positions of signals is correctly inferred, but shrunken
# Less shrinkage as sigma increases, concentrated posteriors with heavy tails
# Issue with beta1 unfortunately. Other signals are recovered. This can be remedied by
# changing the mu vector! Positions and magnitudes can be inferred. 
# Since mu is the mean of a logit, perhaps let it lie bw -3 and 3?
# Lambdas also have a better posterior distribution when mu has structure than when it does not

# Response: use symmetric version and explain what happened to the 1s element

df <- do.call(rbind, results)

df_long <- df %>%
  pivot_longer(
    cols = -sigma_ln,              # Include all columns except sigma_ln
    names_to = "params",           # New column for parameter names
    values_to = "value"            # New column for parameter values
  ) %>%
  mutate(sigma_ln = as.factor(sigma_ln))  # Ensure sigma_ln is a factor

# Plots ----

# beta vs sigma_ln
df_long %>%
  filter(str_detect(params, "^beta") & as.integer(str_extract(params, "\\d+")) <= 9) %>%
  ggplot(aes(y = value, x = params, fill = sigma_ln)) +
  geom_boxplot()+
  #geom_violin(trim = FALSE, draw_quantiles = c(0.025, 0.5, 0.975)) +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Posterior of beta by sigma_ln",
    x = "Sigma_ln Levels",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ sigma_ln, 
             scales = "free_y")  

df %>% 
  ggplot(aes(y = `beta[1]`, x = as.factor(sigma_ln), , fill = as.factor(sigma_ln))) +
  #geom_violin(trim = FALSE, draw_quantiles = c(0.025, 0.5, 0.975) ) +
  geom_boxplot()+
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Distribution of Beta[?] by Sigma_ln",
    x = "Sigma_ln Levels",
    y = ""
  ) +
  #xlim(-5, 5) +   # Restrict y-axis (was x-axis)
  theme_minimal() +
  theme(legend.position = "none")

# logit phi vs sigma_ln
df_long %>%
  filter(str_detect(params, "^logitR2_phi\\[") & 
           as.integer(str_extract(params, "(?<=\\[)\\d+(?=\\])")) <= 9) %>% 
  ggplot(aes(y = value, x = params, fill = sigma_ln)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.025, 0.5, 0.975)) +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Posterior of logit_phi by sigma_ln",
    x = "Sigma_ln Levels",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ sigma_ln, 
             scales = "free_y")  

# lambdas vs sigma_ln
df_long %>%
  filter(str_detect(params, "^lambda\\[") & 
           as.integer(str_extract(params, "(?<=\\[)\\d+(?=\\])")) <= 9) %>% 
  ggplot(aes(y = value, x = params, fill = sigma_ln)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.025, 0.5, 0.975)) +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Posterior of lambda by sigma_ln",
    x = "Sigma_ln Levels",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ sigma_ln, 
             scales = "free_y")  

# R2_phi vs sigma_ln
df_long %>%
  filter(str_detect(params, "^R2_phi\\[") & 
           as.integer(str_extract(params, "(?<=\\[)\\d+(?=\\])")) <= 9) %>% 
  ggplot(aes(y = value, x = params, fill = sigma_ln)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.025, 0.5, 0.975)) +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Posterior of R2_phi by sigma_ln",
    x = "Sigma_ln Levels",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ sigma_ln, 
             scales = "free_y")  # Facet by each beta parameter

# R2_tau2 vs sigma_ln
df_long %>%
  filter( params == "R2_tau2") %>% 
  ggplot(aes(y = value, x = params, fill = sigma_ln)) +
  geom_violin(trim = FALSE, draw_quantiles = c(0.025, 0.5, 0.975)) +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "Posterior of R2_phi by sigma_ln",
    x = "Sigma_ln Levels",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ sigma_ln, 
             scales = "free_y")  # Facet by each beta parameter


# Check the results
df <- do.call(rbind, results)

library(ggplot2)
library(dplyr)

df %>% 
  bayesplot::mcmc_pairs(regex_pars = "^beta\\[[1,2,3,10,15]\\]") 

df %>% 
  bayesplot::mcmc_pairs(regex_pars = "^logitR2_phi\\[[1,2,3,10,15, 19]\\]") 

df %>% 
  bayesplot::mcmc_pairs(regex_pars = "^mu_logitphi\\[[1,2,3,10,15, 19]\\]") 


df %>% 
  bayesplot::mcmc_pairs(regex_pars = "^R2_phi\\[[1,2,3,10,15, 19]\\]") 


df %>% 
  bayesplot::mcmc_areas(regex_pars = "^beta") 

df %>% 
  bayesplot::mcmc_areas(regex_pars = "^mu_logitphi") 

df %>% 
  bayesplot::mcmc_areas(regex_pars = "^logitR2_phi") 

df %>% 
  bayesplot::mcmc_areas(regex_pars = "^R2_phi") 

#------- Relationship between api and sigma_logitphi

# Define the range for x values
x_values <- seq(0.1, 0.5, length.out = 1000)

# Compute the trigamma function for each x value
y_values <- trigamma(x_values)

# Create a data frame for plotting
data <- data.frame(x = x_values, y = y_values)

# Plot using ggplot2
ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Trigamma Function", x = "x", y = "Trigamma(x)") +
  theme_minimal()

trigamma(1e2)
trigamma(50)
trigamma(10)
trigamma(1)
trigamma(0.75) #currently at use
trigamma(0.5) #currently at use
trigamma(0.1) #currently at use 

