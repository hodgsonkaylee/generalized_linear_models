#####################################################
### Simulation function to find regression values ###
#####################################################

sim.data <- function(k,phi,beta,n,nSim){
  # g(mu_t) = sum(x_ti * beta_i) = eta_t, t=1-n, i=1-k
  # randomly draw x values from a uniform distribution with bounds similar to 
  # the coefficients for a linear fit of the dataset
  X1 <- runif(n,0,.50)
  X2 <- runif(n,0,5)
  X3 <- runif(n,0,.60)
  X4 <- runif(n,0,3)
  X <- matrix(c(rep(1,n),X1,X2,X3,X4),nrow=n,ncol=k)
  xbeta <- rep(NA,n)
  mu <- rep(NA,n)
  for(t in 1:n){
    xbeta[t] <- X[t,]%*%beta
    mu[t] <- exp(xbeta[t])/(1+exp(xbeta[t]))
  }
  # find randomly generated y values from different parameter estimates
  p <- matrix(rep(NA,n),nrow=n)
  q <- matrix(rep(NA,n),nrow=n)
  for(t in 1:n) {
    # p and q parameters calculated from mu and phi (functions in article)
    p[t] <- mu[t]*phi
    q[t] <- (1-mu[t])*phi
  }
  # randomly draw the y values from the beta distribution n times and repean nSim times
  # place X in a data frame
  X.df <- data.frame(X)
  # y simulation
  y.sim <- matrix(rep(NA,n*nSim),nrow=n,ncol=nSim)
  # fit the values with both the beta and the linear regression models
  fit <- matrix(rep(NA,nSim*(k+1)),nrow=k+1,ncol=nSim)
  fit.lm <- matrix(rep(NA,nSim*(k+1)),nrow=k,ncol=nSim)
  # Find the mse for both the beta and linear regression models
  mse <- rep(NA,nSim)
  mse.lm <- rep(NA,nSim)
  # Run for-loops to calculate all of the values
  for(l in 1:nSim) {  
    for(t in 1:n){
      y.sim[t,l] <- rbeta(1,p[t],q[t])
    }
  }
  for(l in 1:nSim) {
    fit[,l] <- coef(betareg(y.sim[,l]~X2+X3+X4+X5, data=X.df, link="logit", method="BFGS"))
    fit.lm[,l] <- coef(glm(y.sim[,l]~X2+X3+X4+X5, data=X.df, family=quasibinomial(link = "logit")))
    mse[l] <- mean(betareg(y.sim[,l]~X2+X3+X4+X5, data=X.df, link="logit", method="BFGS")$residuals^2)
    mse.lm[l] <- mean(glm(y.sim[,l]~X2+X3+X4+X5, data=X.df, family=quasibinomial(link = "logit"))$residuals^2)
  }
  y <- apply(y.sim,1,mean)
  fit.dat <- apply(fit,1,mean)
  fit.dat.lm <- apply(fit.lm,1,mean)
  abser <- abs(c(beta,phi) - t(fit.dat))
  abser.lm <- abs(c(beta) - t(fit.dat.lm))
  mseval <- mean(mse)
  mseval.lm <- mean(mse.lm)
  sim.data <- data.frame(matrix(c(X,y,rep(fit.dat,n/6),rep(fit.dat.lm,n/5),
                         rep(abser,n/6),rep(abser.lm,n/5),rep(mse,n),rep(mse.lm,n)),ncol=(k+7),nrow=n))
  return(sim.data)
}

#################################################
### Simulation Function with different values ###
#################################################

# Set number of beta parameters and number of observations for each x 
# (approximately the same as dataset)
k <- 5 
phi <- c(2,5,10,60)
beta <- c(-0.33696,-0.15116,-0.45989,-0.024487,0.46680)
n <- c(150,1200,6000)
nSim <- 1000

# phi=2: X matrix, estimated regression coefficients, absolute error
sim11 <- sim.data(k,phi[1],beta,n[1],nSim)
sim12 <- sim.data(k,phi[1],beta,n[2],nSim)
sim13 <- sim.data(k,phi[1],beta,n[3],nSim)

# phi=5: X matrix, estimated regression coefficients, absolute error
sim21 <- sim.data(k,phi[2],beta,n[1],nSim)
sim22 <- sim.data(k,phi[2],beta,n[2],nSim)
sim23 <- sim.data(k,phi[2],beta,n[3],nSim)

# phi=10: X matrix, estimated regression coefficients, absolute error
sim31 <- sim.data(k,phi[3],beta,n[1],nSim)
sim32 <- sim.data(k,phi[3],beta,n[2],nSim)
sim33 <- sim.data(k,phi[3],beta,n[3],nSim)

# phi=60: X matrix, estimated regression coefficients, absolute error
sim41 <- sim.data(k,phi[4],beta,n[1],nSim)
sim42 <- sim.data(k,phi[4],beta,n[2],nSim)
sim43 <- sim.data(k,phi[4],beta,n[3],nSim)



