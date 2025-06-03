

# Divergence Functions ----------------------------------------------------

# Minimize KL divergence between a ALR Logit normal and a Dirichlet distribution
kl.diri.logitnormal <- function(alpha, ref = 1){
  
  # alpha: concentration vector of the Dirichlet
  
  p <- length(alpha)-1
  mu <- digamma(alpha[-ref])-digamma(alpha[ref])
  #mu <- mu[-ref]
  Sigma <- matrix( trigamma(alpha[ref]),p,p)
  diag(Sigma) <- trigamma(alpha[-ref])+trigamma(alpha[ref])
  
  list(alpha=alpha, mu=mu, Sigma=Sigma)  
  
}

# Stan fit functions ----


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
    return( list(fit.error= TRUE))
  }else{
    return( fit )
  }
  
}

# Generic function to prepare data and call the fitting function
fit_model <- function(model_name, params, additional_params = list()) {
  
  file_model <- file.path("stan/models", paste0(model_name, ".stan"))
  
  # Read the model 
  mod <- cmdstan_model(stan_file = file_model,
                       dir = "stan/compiled", 
                       force_recompile = FALSE)
  
  
  # All models use the same data template. 
  
  dat <- create_data_list(params, additional_params)
  
  run_fit(mod, dat, params)
  
}





