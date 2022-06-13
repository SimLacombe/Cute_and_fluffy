##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##      Perform a Bayesian goodness of 
##      Fit test for the two species
##      occupancy model
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##***************************
## This script contains a set of two functions 
## to test the latent model and the detection model
## Each function works as follows:
## 1. For Nsim, MCMC iterations
## 2. Simulate a dataset with the estimateted parameters
## 3. Get a Chi2 value for the simulated and the true datasets
## 4. Return all chi2 values
##***************************
##

####Latent model-----

get_gof_latent <- function(mm = NULL, data_list = NULL, e = .0001, Nsim = 10000){

  # mm = model_matrix of the model
  # data_list = data list given to the jags file when computing the model
  # e = threshold to calculate the chi2
  # Nsim = number of simulations 
  #-------------------
  # Here the CHI2 is calculated based on the number of times each state is present
  #--------------------
  
  
  #locate each element in data_list to a single variable
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  chi2_l <- chi2.sim_l <- rep(0,Nsim)
   
  cat(paste0("-------------------------------------------------| ", Nsim, "\n"))
  O <- sample(1:nrow(mm),Nsim)
  for(i in 1:Nsim){
    if (!i%%(Nsim/50)) cat("+")
    o <- O[i]
    a <- matrix(mm[o,grep("a\\[", colnames(mm))], ncol = ncov_psi, nrow = nspec)
    b <- matrix(mm[o,grep("b\\[", colnames(mm))], ncol = ncov_gam, nrow = nspec)
    d <- matrix(mm[o,grep("d\\[", colnames(mm))], ncol = ncov_eps, nrow = nspec)
    g <- matrix(mm[o,grep("g\\[", colnames(mm))], ncol = ncov_pi, nrow = nspec)
    h <- matrix(mm[o,grep("h\\[", colnames(mm))], ncol = ncov_tau, nrow = nspec)
    
    z <- matrix(mm[o,grep("z", colnames(mm))], ncol = nseason, nrow = nsite)
      
    psinit  <- invlogit(a %*% t(psi_cov))
    gam     <- invlogit(b %*% t(gam_cov))
    eps     <- invlogit(d %*% t(eps_cov))
    
    gam_one <- invlogit(b %*% t(gam_cov) + g %*% t(pi_cov))
    eps_one <- invlogit(d %*% t(eps_cov) + h %*% t(tau_cov))
    
    # This is the numerator of the softmax function for each of the 4 community
    # states a site can be in.
    fsm <- matrix(0, ncol = 4, nrow = nsite)
    fsm[, 1] <- (1 - psinit[1, ]) * (1 - psinit[2, ]) #------------------------|U
    fsm[, 2] <- psinit[1, ]       * (1 - psinit[2, ]) #------------------------|A
    fsm[, 3] <- (1 - psinit[1, ]) * psinit[2, ]       #------------------------|B
    fsm[, 4] <- psinit[1, ]       * psinit[2, ]       #------------------------|AB

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
    
    #SIMULATED NUMBER OF DETECTIONS : y.sim
    #EXPECTED NUMBER OF DETECTIONS : Y.exp

    z.sim <- array(NA, dim = c(nsite,nseason))
    
    psi <- array(0,dim = c(nsite,nseason,4))
    
    for(m in 1:nsite){
      psi[m,1,] <- fsm[m,]/sum(fsm[m,])
      z.sim[m,1] <- rcat(1,fsm[m,])
      for(t in 2:nseason ){
        psi[m,t,] <- tpm[m,,]%*%psi[m,t-1,]
        psi[m,t,] <- psi[m,t,]/sum(psi[m,t,])
        z.sim[m,t] <- rcat(1,tpm[m, , z.sim[m,t-1]])
      }
    }
    z_freq <- apply(z, 1, function(x) {table(factor(x, levels = 1:4))/length(which(!is.na(x)))})
    z.sim_freq <- apply(z.sim, 1, function(x) {table(factor(x, levels = 1:4))/length(which(!is.na(x)))})
    z.exp_freq <- apply(psi, 1, function(x) {apply(x,2,function(x) sum(x, na.rm=T))/sum(x, na.rm=T)})
    
  #COMPUTE CHI²
    
    chi2 <- sum((z_freq - z.exp_freq)^2/(z.exp_freq+e))
    chi2.sim <- sum((z.sim_freq - z.exp_freq)^2/(z.exp_freq+e))
    chi2_l[i] <- chi2
    chi2.sim_l[i] <- chi2.sim

  }
  return(data.frame(chi2.obs = chi2_l, chi2.sim=chi2.sim_l))
}

####Observation model-----

get_gof_observation <- function(mm = NULL, data_list = NULL, e = .0001, Nsim = 10000){
  
  # mm = model_matrix of the model
  # data_list = data list given to the jags file when computing the model
  # e = threshold to calculate the chi2
  # Nsim = number of simulations 
  #-------------------
  # Here the CHI2 is calculated based on the number of time each species is observed
  # We also return a data frame where the chi2 are grouped by species, site, or covariate
  #--------------------
  
  
  #locate each element in data_list to a single variable
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  detfreq <- apply(y, c(1,2), function(x) {table(factor(x, levels = 1:4))})
  
  chi2_l <- chi2.sim_l <- rep(0,Nsim)
  chi2_RF_l <- chi2_RF.sim_l <- chi2_AF_l <- chi2_AF.sim_l <- rep(0,Nsim)
  chi2_sites_m <- chi2_sites.sim_m <- matrix(0, nsite, Nsim)  
  
  cat(paste0("-------------------------------------------------| ", Nsim, "\n"))
  O <- sample(1:nrow(mm),Nsim)
  for(i in 1:Nsim){
    if (!i%%(Nsim/50)) cat("+")
    o <- O[i]
    
    f <- matrix(mm[o,grep("f\\[", colnames(mm))], ncol = ncov_rho, nrow = nspec)
    yr_rho <- matrix(mm[o, grep("yr_rho", colnames(mm))], ncol = nyear, nrow = nspec)
    bait_ <- matrix(mm[o, grep("bait_", colnames(mm))], ncol = nsite, nrow = nspec)
    
    # estimated community state at time t and site j
    z <- matrix(mm[o,grep("z", colnames(mm))], ncol = nseason, nrow = nsite)
    
    rho       <- invlogit(yr_rho[,year_cov] + f %*% t(rho_cov))
    rho_bait  <- invlogit(yr_rho[,year_cov] + f %*% t(rho_cov) + bait_)
    
    rdm <- array(0, dim = c(nsite, 4, 4,2))
    ## Bait absent
    # TS = U
    rdm[, 1, 1, 1] <- 1 #--------------------------------------------------|OS = U
    rdm[, 2, 1, 1] <- 0 #--------------------------------------------------|OS = A
    rdm[, 3, 1, 1] <- 0 #--------------------------------------------------|OS = B
    rdm[, 4, 1, 1] <- 0 #--------------------------------------------------|OS = AB
    # TS = A
    rdm[, 1, 2, 1] <- (1 - rho[1, ])#--------------------------------------|OS = U
    rdm[, 2, 2, 1] <-      rho[1, ] #--------------------------------------|OS = A
    rdm[, 3, 2, 1] <- 0 #--------------------------------------------------|OS = B
    rdm[, 4, 2, 1] <- 0 #--------------------------------------------------|OS = AB
    # TS = B
    rdm[, 1, 3, 1] <- (1 - rho[2, ]) #-------------------------------------|OS = U
    rdm[, 2, 3, 1] <- 0 #--------------------------------------------------|OS = A
    rdm[, 3, 3, 1] <-      rho[2, ]  #-------------------------------------|OS = B
    rdm[, 4, 3, 1] <- 0 #--------------------------------------------------|OS = AB
    # TS = AB
    rdm[, 1, 4, 1] <- (1 - rho[1, ]) * (1 - rho[2, ]) #-------------------|OS = U
    rdm[, 2, 4, 1] <-      rho[1, ]  * (1 - rho[2, ]) #-------------------|OS = A
    rdm[, 3, 4, 1] <- (1 - rho[1, ]) *      rho[2, ]  #-------------------|OS = B
    rdm[, 4, 4, 1] <-      rho[1, ]  *      rho[2, ]  #-------------------|OS = AB
    
    ## Bait present
    # TS = U
    rdm[, 1, 1, 2] <- 1 #--------------------------------------------------|OS = U
    rdm[, 2, 1, 2] <- 0 #--------------------------------------------------|OS = A
    rdm[, 3, 1, 2] <- 0 #--------------------------------------------------|OS = B
    rdm[, 4, 1, 2] <- 0 #--------------------------------------------------|OS = AB
    # TS = A
    rdm[, 1, 2, 2] <- (1 - rho_bait[1, ])#---------------------------------|OS = U
    rdm[, 2, 2, 2] <-      rho_bait[1, ] #---------------------------------|OS = A
    rdm[, 3, 2, 2] <- 0 #--------------------------------------------------|OS = B
    rdm[, 4, 2, 2] <- 0 #--------------------------------------------------|OS = AB
    # TS = B
    rdm[, 1, 3, 2] <- (1 - rho_bait[2, ]) #--------------------------------|OS = U
    rdm[, 2, 3, 2] <- 0 #--------------------------------------------------|OS = A
    rdm[, 3, 3, 2] <-      rho_bait[2, ]  #--------------------------------|OS = B
    rdm[, 4, 3, 2] <- 0 #--------------------------------------------------|OS = AB
    # TS = AB
    rdm[, 1, 4, 2] <- (1 - rho_bait[1, ]) * (1 - rho_bait[2, ]) #----------|OS = U
    rdm[, 2, 4, 2] <-      rho_bait[1, ]  * (1 - rho_bait[2, ]) #----------|OS = A
    rdm[, 3, 4, 2] <- (1 - rho_bait[1, ]) *      rho_bait[2, ]  #----------|OS = B
    rdm[, 4, 4, 2] <-      rho_bait[1, ]  *      rho_bait[2, ]  #----------|OS = AB
    
    #SIMULATED NUMBER OF DETECTIONS : y.sim
    
    y.sim <- array(NA, dim = c(nsite,nseason, nsurvey))
    y.exp <- array(NA, dim = c(nsite,nseason, nsurvey, 4))
    for(m in 1:nsite){
      for(t in 1:nseason ){
        for(k in 1:nsurvey){
          if(!is.na(y[m,t,k])){
            bait_ <- bait[m,t,k]
            y.sim[m,t,k] <- rcat(1, rdm[m, ,z[m,t], bait_])
            y.exp[m,t,k,] <- rdm[m, ,z[m,t], bait_]/sum(rdm[m, ,z[m,t], bait_])
          }
        }
      }
    }

    detfreq_sim <- apply(y.sim, c(1,2), function(x) {table(factor(x, levels = 1:4))})
    detfreq_exp <- apply(y.exp, c(1,2), function(x) {apply(x,2,function(x) sum(x, na.rm=T))})
    
    #Get chi 2 residuals 
    
    chi2_eps <- (detfreq - detfreq_exp)^2/(detfreq_exp+e)
    chi2.sim_eps <- (detfreq_sim - detfreq_exp)^2/(detfreq_exp+e)
    
    #COMPUTE global CHI²
    
    chi2 <- sum(chi2_eps)
    chi2.sim <- sum(chi2.sim_eps)
    chi2_l[i] <- chi2
    chi2.sim_l[i] <- chi2.sim
    
    #COMPUTE CHI² for species 
    
    chi2_RF <- sum(chi2_eps[c(2,4),,])
    chi2_RF.sim <- sum(chi2.sim_eps[c(2,4),,])
    
    chi2_AF <- sum(chi2_eps[c(3,4),,])
    chi2_AF.sim <- sum(chi2.sim_eps[c(3,4),,])
    
    chi2_RF_l[i] <- chi2_RF
    chi2_RF.sim_l[i] <- chi2_RF.sim
    
    chi2_AF_l[i] <- chi2_AF
    chi2_AF.sim_l[i] <- chi2_AF.sim
    
    #COMPUTE CHI² for sites 
    
    chi2_sites <- apply(chi2_eps, 2, sum)
    chi2_sites.sim <- apply(chi2.sim_eps, 2, sum)
    
    chi2_sites_m[, i] <- chi2_sites
    chi2_sites.sim_m[, i] <- chi2_sites.sim
  }
  
  res_sites_df <- as.data.frame(t(chi2_sites_m - chi2_sites.sim_m))
  colnames(res_sites_df) <- covs$site.year
  
  res_sites_df <- res_sites_df%>%
    pivot_longer(cols = covs$site.year)%>%
    mutate(year      = rep(covs$year,           1000),
           loc       = rep(covs$loc,            1000),
           treatment = rep(covs$treatment,      1000),
           region    = rep(substr(covs$loc,1,1), 1000),
           CLG       = rep(covs$CLG,            1000),
           FTG       = rep(covs$FTG,            1000),
           rodents   = rep(covs$rodents_fall,   1000))%>%
    mutate(CLG.bins      = cut_number(CLG,     n= 10),
           FTG.bins      = cut_number(FTG,     n= 10),
           rodents.bins  = cut_number(rodents, n= 10))
  
  return(list(GLOBAL = data.frame(chi2.obs = chi2_l,
                                  chi2.sim=chi2.sim_l),
              SPECIES = data.frame(chi2.obs = c(chi2_RF_l,chi2_AF_l),
                                   chi2.sim = c(chi2_RF.sim_l,chi2_AF.sim_l), 
                                   species = rep(c("RF","AF"), each = Nsim)),
              SITES = res_sites_df))
}

