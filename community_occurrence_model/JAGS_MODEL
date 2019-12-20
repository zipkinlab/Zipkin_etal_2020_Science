model {
  #Define prior distributions for community-level model parameters
  for (t in 1:nPeriod){
    omega[t] ~ dunif(0,1)
    u.mu[t] ~ dnorm(0, 0.368)
  }#t
  v.mu ~ dnorm(0, 0.368)  
  u.sigma ~ dgamma(0.1, 0.1)
  u.tau <- 1/(u.sigma*u.sigma)
  rho ~ dnorm(0,0.368)
  
  #Specify the likelihood
  
  for (i in 1:nSum){
    for (t in 1:nPeriod){
      
      #Define the species-level priors
      w[i,t] ~ dbern(omega[t])  
      u[i,t] ~ dnorm(u.mu[t], u.tau)
     
      #Create a loop to estimate the Z matrix (true occurrence for species i)
      #at point j.
      for (j in 1:nTran){
        logit(psi[i,j,t]) <- u[i,t] 
        psi.mu[i,j,t] <- psi[i,j,t]*w[i,t]
        Z[i,j,t] ~ dbern(psi.mu[i,j,t])
      }#j
    }#t

    for (k in 1:nSamples) {
      logit(p[k,i]) <-  v.mu + rho*(u[i,period[k]]-u.mu[period[k]])
      mu.p[k,i] <- p[k,i]*Z[i,transects[k],period[k]]
      X[k,i] ~ dbern(mu.p[k,i])
      
      #Create simulated dataset to calculate the Bayesian p-value
      Xnew[k,i] ~ dbern(mu.p[k,i])
      d[k,i]<-  abs(X[k,i] - mu.p[k,i]) 
      dnew[k,i]<- abs(Xnew[k,i] - mu.p[k,i]) 
      d2[k,i]<- pow(d[k,i],2)  
      dnew2[k,i]<- pow(dnew[k,i],2) 
    }#k
  }#i
  
  #Calculate the discrepancy measure, defined as the mean(p.fit > p.fitnew) 
  #using only observed species (1 through 36)
  p.fit<-sum(d2[1:nSamples,1:36]) 
  p.fitnew<-sum(dnew2[1:nSamples,1:36])
  p.diff<-step(p.fit - p.fitnew)

}
