##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##      get stationary occupancy probabilities
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##

library(doParallel)

get_steady_state <- function(mm = NULL, data_list = NULL, param_covs, cov, Nsim = 1000, ncores = 4){
  
  mm <- mm[, -grep("z", colnames(mm))]
  
  covs.m <- cbind(1,0,0,0,rep(c(0,1),each=100))
  colnames(covs.m) <- c("int", "CLG", "FTG", "rodents_fall", "treatment")
  covs.m[, cov[1]] <- rep(seq(-4,4,length.out = 100),2)
    
  #locate each element in data_list to a single variable
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  nsite <- nrow(covs.m)
   
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  n <- nrow(mm)

  O <- sample(1:nrow(mm),Nsim)
  
  steady.array <- foreach (i = 1:Nsim,
                           .packages = c('LaplacesDemon', 'expm'))%dopar%{
    o <- O[i]
    
    steady.m <- array(0, dim = c(nsite, 4))
    
    a <- matrix(mm[o,grep("a\\[", colnames(mm))], ncol = ncov_psi, nrow = nspec)
    b <- matrix(mm[o,grep("b\\[", colnames(mm))], ncol = ncov_gam, nrow = nspec)
    d <- matrix(mm[o,grep("d\\[", colnames(mm))], ncol = ncov_eps, nrow = nspec)
    g <- matrix(mm[o,grep("g\\[", colnames(mm))], ncol = ncov_pi, nrow = nspec)
    h <- matrix(mm[o,grep("h\\[", colnames(mm))], ncol = ncov_tau, nrow = nspec)
    
      
    psinit  <- invlogit(a %*% t(covs.m[,param_covs[["psi"]]]))
    gam     <- invlogit(b %*% t(covs.m[,param_covs[["gam"]]]))
    eps     <- invlogit(d %*% t(covs.m[,param_covs[["eps"]]]))
    
    gam_one <- invlogit(b %*% t(covs.m[,param_covs[["gam"]]]) + g %*% t(covs.m[,param_covs[["pi"]]]))
    eps_one <- invlogit(d %*% t(covs.m[,param_covs[["eps"]]]) + h %*% t(covs.m[,param_covs[["tau"]]]))
    
    # Transition probabilities for each site
    tpm <- array(0, dim = c(nsite, 4, 4))
    # U to ...
    tpm[, 1, 1] <- (1 - gam[1, ])     * (1 - gam[2, ])     #-----------------|U
    tpm[, 2, 1] <- gam[1, ]           * (1 - gam_one[2, ]) #-----------------|A
    tpm[, 3, 1] <- (1 - gam[1, ])     * gam[2, ]           #-----------------|B
    tpm[, 4, 1] <- gam[1,  ]          * gam_one[2, ]       #-----------------|AB
    # A to ...
    tpm[, 1, 2] <- eps[1, ]           * (1 - gam_one[2, ]) #-----------------|U
    tpm[, 2, 2] <- (1 - eps[1, ])     * (1 - gam_one[2, ]) #-----------------|A
    tpm[, 3, 2] <- eps_one[1, ]       * gam_one[2, ]       #-----------------|B
    tpm[, 4, 2] <- (1 - eps_one[1, ]) * gam_one[2, ]       #-----------------|AB 
    # B to ...
    tpm[, 1, 3] <- (1 - gam_one[1, ]) * eps[2, ]           #-----------------|U
    tpm[, 2, 3] <- gam_one[1, ]       * eps_one[2, ]       #-----------------|A
    tpm[, 3, 3] <- (1 - gam_one[1, ]) * (1 - eps[2, ])     #-----------------|B
    tpm[, 4, 3] <- gam_one[1, ]       * (1 - eps_one[2, ]) #-----------------|AB
    # AB to ..
    tpm[, 1, 4] <- eps_one[1, ]       * eps_one[2, ]       #-----------------|U
    tpm[, 3, 4] <- (1 - eps_one[1, ]) * eps_one[2, ]       #-----------------|A
    tpm[, 2, 4] <- eps_one[1, ]       * (1 - eps_one[2, ]) #-----------------|B
    tpm[, 4, 4] <- (1 - eps_one[1, ]) * (1 - eps[2, ])     #-----------------|AB

    for(i in 1:nsite){
      mark <- new("markovchain", states = c("U", "A", "B", "AB"),byrow = TRUE, t(tpm[i,,]))
      steady.m[i,]  <-   steadyStates(mark)
    }
    steady.m
  }%>%sapply(.,I,simplify="array")
  stopCluster(cl)
  return(steady.array)
}