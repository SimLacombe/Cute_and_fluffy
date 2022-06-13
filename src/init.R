##############################
#
#   Get initial parametter values 
#   to run the occupancy model
#
##############################

# generates initial values for mcmc chains for all parameters

inits <- function(chain){
  gen_list <- function(chain = chain){
    l.init <- list(z = init,
                   a = matrix(rnorm(data_list$nspec * data_list$ncov_psi),
                             ncol = data_list$ncov_psi, nrow = data_list$nspec),
                   b = matrix(rnorm(data_list$nspec * data_list$ncov_gam),
                              ncol = data_list$ncov_gam, nrow = data_list$nspec),
                   d = matrix(rnorm(data_list$nspec * data_list$ncov_eps),
                             ncol = data_list$ncov_eps, nrow = data_list$nspec),
                   f = matrix(rnorm(data_list$nspec * data_list$ncov_rho),
                             ncol = data_list$ncov_rho, nrow = data_list$nspec),
                  g = matrix(rnorm(data_list$nspec * data_list$ncov_pi),
                             ncol = data_list$ncov_pi, nrow = data_list$nspec),
                  h = matrix(rnorm(data_list$nspec * data_list$ncov_tau),
                             ncol = data_list$ncov_tau, nrow = data_list$nspec),
                  yr_rho = matrix(rnorm(data_list$nspec*data_list$nyear),
                             nrow = data_list$nspec, ncol = data_list$nyear),
                  bait_ = rnorm(data_list$nspec),
                  .RNG.name = switch(chain,
                                    "1" = "base::Wichmann-Hill",
                                    "2" = "base::Marsaglia-Multicarry",
                                    "3" = "base::Super-Duper",
                                    "4" = "base::Mersenne-Twister",
                                    "5" = "base::Wichmann-Hill",
                                    "6" = "base::Marsaglia-Multicarry",
                                    "7" = "base::Super-Duper",
                                    "8" = "base::Mersenne-Twister"),
                    .RNG.seed = sample(1:1e+06, 1))
    return(l.init)
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}