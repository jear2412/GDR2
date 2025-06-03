# Helper function to create list configurations
create_fit_params <- function(id, mcmc_name, fit, extra_params) {
  list(
    id = id,
    mcmc_name = mcmc_name,
    fit = fit,
    extra_params = extra_params
  )
}

r2d2_params <- function(mean_R2, prec_R2, api) {
  list(R2_mean = mean_R2, R2_prec = prec_R2, gen_R2_alpha_function = "R2_alpha_api", api = api)
}

#--- Models 

betaprime_list <- create_fit_params( "BP", "betaprime", "betaprime", list(a = 0.5, b = 1))

dirlaplace_list <- create_fit_params("DL", "dirichlet_laplace", "dirichletlaplace", list(api = 0.5))

horseshoe_list <- create_fit_params("HS", "horseshoe", "horseshoe", list(NULL))

normalgamma_list <- create_fit_params("NG", "normalgamma", "normalgamma", list(NULL))

logitr2d2_params <- function(mean_R2, prec_R2, api,  inference_config) {
  list(R2_mean = mean_R2, 
       R2_prec = prec_R2,
       gen_R2_alpha_function = "R2_alpha_api",
       api = api, 
       inference_config = inference_config)
}

# R2D2 models

r2d2_listdd <- create_fit_params("D2dd", "r2d2", "r2d2", r2d2_params(0.5, 1, 0.5))
r2d2_listud <- create_fit_params("D2ud", "r2d2", "r2d2", r2d2_params(0.5, 2, 0.5))
r2d2_listdu <- create_fit_params("D2du", "r2d2", "r2d2", r2d2_params(0.5, 1, 1))
r2d2_listuu <- create_fit_params("D2uu", "r2d2", "r2d2", r2d2_params(0.5, 2, 1))


# logit R2 models

# logit R2 KL div

logitr2d2_list_dd <- create_fit_params("L2dd", "logitr2", "logitr2", 
                                     logitr2d2_params(0.5, 1, 0.5,
                                                      inference_config = list( method = "kl_divergence" )))

logitr2d2_list_ud <- create_fit_params("L2ud", "logitr2", "logitr2", 
                                     logitr2d2_params(0.5, 2, 0.5,
                                                      inference_config = list( method = "kl_divergence" )))


logitr2d2_list_du <- create_fit_params("L2du", "logitr2", "logitr2", 
                                         logitr2d2_params(0.5, 1, 1,
                                                          inference_config = list( method = "kl_divergence" )))

logitr2d2_list_uu <- create_fit_params("L2uu", "logitr2", "logitr2", 
                                         logitr2d2_params(0.5, 2, 1,
                                                          inference_config = list( method = "kl_divergence" )))

# logit R2 KL scales 

logitr2d2_listd_dd_s <- create_fit_params("L2dd_s", "logitr2", "logitr2", 
                                     logitr2d2_params(0.5, 1, 0.5,
                                                    inference_config = list( method = "kl_divergence_scales" )))

logitr2d2_listd_ud_s <- create_fit_params("L2ud_s", "logitr2", "logitr2", 
                                     logitr2d2_params(0.5, 2, 0.5,
                                                    inference_config = list( method = "kl_divergence_scales" )))

logitr2d2_listd_du_s <- create_fit_params("L2du_s", "logitr2", "logitr2", 
                                           logitr2d2_params(0.5, 1, 1,
                                                            inference_config = list( method = "kl_divergence_scales" )))

logitr2d2_listd_uu_s <- create_fit_params("L2uu_s", "logitr2", "logitr2", 
                                           logitr2d2_params(0.5, 2, 1,
                                                            inference_config = list( method = "kl_divergence_scales" )))
