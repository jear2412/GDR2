
# Helper function to create the common data list
create_data_list <- function(params, additional_params = list()) {
  N <- params$n
  p <- params$p
  X <- params$X
  y <- params$y
  Ntest <- params$ntest
  Xtest <- params$Xtest
  ytest <- params$ytest
  scale_sigma <- sd(y)
  
  # Base data list
  dat <- list(
    N = N,
    p = p + 1,
    X = cbind(rep(1, N), X),
    y = as.numeric(y),
    Ntest = Ntest,
    Xtest = cbind(rep(1, Ntest), Xtest),
    ytest = as.numeric(ytest),
    scale_sigma = scale_sigma,
    prior_only = params$prior_only
  )
  
  # Add additional parameters
  dat <- modifyList(dat, additional_params)
  
  return(dat)
}

# Helper function to perform the sampling and error handling
run_fit <- function(mod, dat, params){
  
  fit <- try(
    mod$sample(
      data = dat,
      seed = params$seed,
      refresh = 500,
      chains = params$nchains,
      iter_warmup = params$iter_warmup,
      iter_sampling = params$iter_sample,
      adapt_delta = params$adapt_delta 
      #init = params$init , 
    ),
    silent = FALSE
  )
  
  # TODO: Improve this
  if (is(fit, "try-error")) {
    fit.error <- TRUE
    fit <- NULL
  } else {
    fit.error <- fit$return_codes()  
  }
  
  if(sum(fit.error) == length(fit.error)) {
    #In case there is an error
    return(list(fit.error= TRUE))
  }else{
    return(list(fit= fit, 
                fit.error= FALSE))
  }
  
}

# Generic function to prepare data and call the fitting function
fit_model <- function(model_path, params, additional_params = list()) {
  
  file_model <- file.path("stan/models", paste0(model_path, ".stan"))
  
  # Read the model 
  mod <- cmdstan_model(stan_file = file_model,
                       dir = "stan/compiled", 
                       force_recompile = FALSE)
  
  
  # All models use the same data template.You might add new data that
  # is specific to the fit.
  
  dat <- create_data_list(params, additional_params)
  
  run_fit(mod, dat, params)
  
}


#-- Fit specific information ----

betaprimefit <- function(params){
  mod <- "betaprime"
  additional_params <- list(a = params$a, b = params$b)
  fit_model(mod, params, additional_params)
}



#-- Dirichlet Laplace

dirichletlaplacefit <- function(params){
  mod <-  "dirichletlaplace"
  additional_params <- list(alpha = params$alpha)
  fit_model(mod, params, additional_params)
  
}

#-- logitr2d2
logitr2fit <- function(params){
  mod_dir <- "" # will be modified below
  flag <- params$method
  
  # Extract common parameters
  R2_mean <- params$R2_mean
  R2_prec <- params$R2_prec
  
  # Determine model and additional parameters based on method
  additional_params <- switch(flag,
                              "kl_divergence" = {
                                list(mod = paste0(mod_dir, "logitR2"), 
                                     additional = list(mu_logitphi = params$mu_logitphi, 
                                                       chol_sigma_logitphi = params$chol_sigma_logitphi))
                              },
                              "kl_divergence_scales" = {
                                list(mod = paste0(mod_dir, "logitR2"), 
                                     additional = list(mu_logitphi = params$mu_logitphi, 
                                                       chol_sigma_logitphi = params$chol_sigma_logitphi))
                              },
                              "user_specified" = {
                                list(mod = paste0(mod_dir, "logitR2"), 
                                     additional = list(mu_logitphi = params$mu_logitphi, 
                                                       chol_sigma_logitphi = params$chol_sigma_logitphi))
                              },
                              stop("Invalid inference method specified in `params$method`")
  )
  
  
  # Fit the model
  fit_model(
    model_path = additional_params$mod,
    params = params,
    additional_params = c(additional_params$additional, 
                          list(R2_mean = R2_mean, 
                               R2_prec = R2_prec))
  )
}



#-- Horseshoe
horseshoefit <- function(params){
  
  mod <-  "horseshoe"
  
  additional_params <- list()
  
  fit_model(mod, params, additional_params)
  
}


#-- Normal Gamma
normalgammafit <- function(params){
  mod <-  "normalgamma"
  
  p <- params$p
  M <- 2*0.5*sqrt(log(2))*p #see paper 
  
  additional_params <- list( M = M)
  
  fit_model(mod, params, additional_params)
  
}


#-- R2D2
r2d2fit <- function(params){
  mod <-  "R2D2"
  
  #fit-specific params
  R2_mean =  params$R2_mean
  R2_prec =  params$R2_prec
  R2_alpha = params$R2_alpha
  
  additional_params <- list(R2_mean = R2_mean, 
                            R2_prec = R2_prec, 
                            R2_alpha = R2_alpha)
  
  fit_model(mod, params, additional_params)
  
  
}


#-- Regularized Horseshoe
rhorseshoefit <- function(params){
  
  mod <-  "rhorseshoe"
  
  
  hs_df =  params$hs_df
  hs_df_global = params$hs_df_global
  hs_df_slab = params$hs_df_slab
  hs_scale_slab = params$hs_scale_slab
  p0 = params$p0
  
  additional_params <- list( hs_df = hs_df, 
                             hs_df_global = hs_df_global,
                             hs_df_slab = hs_df_slab, 
                             hs_scale_slab = hs_scale_slab, 
                             p0 = p0)
  
  fit_model(mod, params, additional_params)
  
  
}

