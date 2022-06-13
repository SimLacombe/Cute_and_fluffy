
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##      JAGS file to run the model
##      called in fit_model.R
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##


model{
  ######
  # Detection model
  ######
  for(j in 1:nsite) {
    for(ti in 1:nseason) {
      for(day in 1:nsurvey) {
        y[j, ti, day] ~ dcat( rdm[j, (1:nout) , z[j, ti], bait[j,ti, day] ] )
      }
    }
  }
  #######
  # Latent state model
  #######
  for(j in 1:nsite) {
    # for first season
    z[j, 1] ~ dcat( fsm[j, ( 1:nout )] )
    for(t in 2:nseason){
      z[j, t] ~ dcat( tpm[j, ( 1:nout ) , z[ j, t-1]] )
    }
  }
  for( j in 1:nsite ) {
    ######
    # Fill in all of the transition probabilities
    ######
    ######
    # Latent state for first season, fsm = transition vector for season 1
    ######
    # First season probabilities for each state
    fsm[j, 1] <- (1 - psinit[1, j]) * (1 - psinit[2, j]) #------------------------|U
    fsm[j, 2] <- psinit[1, j]       * (1 - psinit[2, j]) #------------------------|A
    fsm[j, 3] <- (1 - psinit[1, j]) * psinit[2, j]       #------------------------|B
    fsm[j, 4] <- psinit[1, j]       * psinit[2, j]       #------------------------|AB
    ######
    # Latent state for dynamic part of model
    # tpm = transition probability matrix. All columns sum to 1.
    # dim(tpm)[1] = site j 
    # dim(tpm)[2] = state at time t
    # dim(tpm)[3] = state at time t-1
    ######
    # U to ...
    tpm[j, 1, 1] <- (1 - gam[1, j])     * (1 - gam[2, j])     #-----------------|U
    tpm[j, 2, 1] <- gam[1, j]           * (1 - gam_one[2, j]) #-----------------|A
    tpm[j, 3, 1] <- (1 - gam[1, j])     * gam[2, j]           #-----------------|B
    tpm[j, 4, 1] <- gam[1,  j]          * gam_one[2, j]       #-----------------|AB
    # A to ...
    tpm[j, 1, 2] <- eps[1, j]           * (1 - gam_one[2, j]) #-----------------|U
    tpm[j, 2, 2] <- (1 - eps[1, j])     * (1 - gam_one[2, j]) #-----------------|A
    tpm[j, 3, 2] <- eps_one[1, j]       * gam_one[2, j]       #-----------------|B
    tpm[j, 4, 2] <- (1 - eps_one[1, j]) * gam_one[2, j]       #-----------------|AB 
    # B to ...
    tpm[j, 1, 3] <- (1 - gam_one[1, j]) * eps[2, j]           #-----------------|U
    tpm[j, 2, 3] <- gam_one[1, j]       * eps_one[2, j]       #-----------------|A
    tpm[j, 3, 3] <- (1 - gam_one[1, j]) * (1 - eps[2, j])     #-----------------|B
    tpm[j, 4, 3] <- gam_one[1, j]       * (1 - eps_one[2, j]) #-----------------|AB
    # AB to ..
    tpm[j, 1, 4] <- eps_one[1, j]       * eps_one[2, j]       #-----------------|U
    tpm[j, 3, 4] <- (1 - eps_one[1, j]) * eps_one[2, j]       #-----------------|A
    tpm[j, 2, 4] <- eps_one[1, j]       * (1 - eps_one[2, j]) #-----------------|B
    tpm[j, 4, 4] <- (1 - eps_one[1, j]) * (1 - eps[2, j])     #-----------------|AB
    
    ######
    # detection matrix (OS = observed state, TS = true state)
    # rdm = rho detection matrix. Each row sums to 1.
    # OS along rows, TS along columns
    ######
    ## Bait absent
    # TS = U
    rdm[j, 1, 1, 1] <- 1 #----------------------------------------------------|OS = U
    rdm[j, 2, 1, 1] <- 0 #----------------------------------------------------|OS = A
    rdm[j, 3, 1, 1] <- 0 #----------------------------------------------------|OS = B
    rdm[j, 4, 1, 1] <- 0 #----------------------------------------------------|OS = AB
    # TS = A
    rdm[j, 1, 2, 1] <- (1 - rho[1, j])#---------------------------------------|OS = U
    rdm[j, 2, 2, 1] <-      rho[1, j] #---------------------------------------|OS = A
    rdm[j, 3, 2, 1] <- 0 #----------------------------------------------------|OS = B
    rdm[j, 4, 2, 1] <- 0 #----------------------------------------------------|OS = AB
    # TS = B
    rdm[j, 1, 3, 1] <- (1 - rho[2, j]) #--------------------------------------|OS = U
    rdm[j, 2, 3, 1] <- 0 #----------------------------------------------------|OS = A
    rdm[j, 3, 3, 1] <-      rho[2, j]  #--------------------------------------|OS = B
    rdm[j, 4, 3, 1] <- 0 #----------------------------------------------------|OS = AB
    # TS = AB
    rdm[j, 1, 4, 1] <- (1 - rho[1, j]) * (1 - rho[2, j]) #--------------------|OS = U
    rdm[j, 2, 4, 1] <-      rho[1, j]  * (1 - rho[2, j]) #--------------------|OS = A
    rdm[j, 3, 4, 1] <- (1 - rho[1, j]) *      rho[2, j]  #--------------------|OS = B
    rdm[j, 4, 4, 1] <-      rho[1, j]  *      rho[2, j]  #--------------------|OS = AB
    
    ## Bait present
    # TS = U
    rdm[j, 1, 1, 2] <- 1 #----------------------------------------------------|OS = U
    rdm[j, 2, 1, 2] <- 0 #----------------------------------------------------|OS = A
    rdm[j, 3, 1, 2] <- 0 #----------------------------------------------------|OS = B
    rdm[j, 4, 1, 2] <- 0 #----------------------------------------------------|OS = AB
    # TS = A
    rdm[j, 1, 2, 2] <- (1 - rho_bait[1, j])#----------------------------------|OS = U
    rdm[j, 2, 2, 2] <-      rho_bait[1, j] #----------------------------------|OS = A
    rdm[j, 3, 2, 2] <- 0 #----------------------------------------------------|OS = B
    rdm[j, 4, 2, 2] <- 0 #----------------------------------------------------|OS = AB
    # TS = B
    rdm[j, 1, 3, 2] <- (1 - rho_bait[2, j]) #---------------------------------|OS = U
    rdm[j, 2, 3, 2] <- 0 #----------------------------------------------------|OS = A
    rdm[j, 3, 3, 2] <-      rho_bait[2, j]  #---------------------------------|OS = B
    rdm[j, 4, 3, 2] <- 0 #----------------------------------------------------|OS = AB
    # TS = AB
    rdm[j, 1, 4, 2] <- (1 - rho_bait[1, j]) * (1 - rho_bait[2, j]) #----------|OS = U
    rdm[j, 2, 4, 2] <-      rho_bait[1, j]  * (1 - rho_bait[2, j]) #----------|OS = A
    rdm[j, 3, 4, 2] <- (1 - rho_bait[1, j]) *      rho_bait[2, j]  #----------|OS = B
    rdm[j, 4, 4, 2] <-      rho_bait[1, j]  *      rho_bait[2, j]  #----------|OS = AB
    
    ######
    # Fill in the linear predictors for the transition matrices
    ######
    for( i in 1:nspec ) {
      # base occupancy
      psinit[i ,j]  <- ilogit(inprod( a[i, ], psi_cov[j, ] ))
      # base colonization
      gam[i, j] <-     ilogit(inprod( b[i, ], gam_cov[j, ] ))
      # base extinction
      eps[i, j] <-     ilogit(inprod( d[i, ], eps_cov[j, ] ))
      # base detection probability
      rho[i,j] <-      ilogit(yr_rho[i, year_cov[j,1]] + inprod(f[i, ], rho_cov[j, ] )) 
      rho_bait[i,j] <- ilogit(yr_rho[i, year_cov[j,1]] + inprod(f[i, ], rho_cov[j, ] ) + bait_[i]) 
      # inxs on colonization
      pi[i,j] <- inprod( g[i, ], pi_cov[j, ] )
      #inxs on extinction
      tau[i,j] <- inprod( h[i, ], pi_cov[j, ] )
      # linear predictor for colonization | one species present at t-1
      gam_one[i, j] <- ilogit(inprod( b[i, ], gam_cov[j, ] ) + pi[i, j])
      # linear predictor for extinction | one species present at t-1
      eps_one[i, j] <- ilogit(inprod( d[i, ], eps_cov[j, ] ) + tau[i, j])
    }
  } # closes for loop for j (sites) all the way up at the top of the model
  
  #####
  # Priors
  ######
  for(i in 1:nspec){
    # Initial Occupancy
    for( psip in 1:ncov_psi ){
      a[i, psip] ~ dlogis(0, 1)
    }
    # Colonization
    for( gamp in 1:ncov_gam ){
      b[i, gamp] ~ dlogis(0, 1)
    }  
    # Extinction
    for( epsp in 1:ncov_eps ){
      d[i, epsp] ~ dlogis(0, 1)
    }
    # Detection
    for( rhop in 1:ncov_rho ){
      f[i, rhop] ~ dlogis(0, 1)
    }
    # Inxs on colonization
    for( pip in 1:ncov_pi ){
      g[i, pip] ~ dlogis(0, 1)
    }
    # Inxs on extinction
    for( taup in 1:ncov_tau ){
      h[i, taup] ~ dlogis(0, 1)
    }
    for ( yrp in 1:nyear ){
      yr_rho[i, yrp] ~ dlogis(0, 1)
    }
    bait_[i] ~ dlogis(0,1)
  }
}