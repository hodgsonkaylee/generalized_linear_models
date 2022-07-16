mu <- 5
beta <- 2
theta <- c(mu,beta)
n <- 100
B <- 1000

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x <- replicate(B, Finv(runif(n),mu,beta))
x <- matrix(x, nrow=B, ncol=n, byrow=TRUE)

# MLE
mle <- function(x){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle <- mle(x[1,])
theta.mle <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle[i,] <- mle(x[i,])
}
theta.mle.mean <- apply(theta.mle,2,mean)
bias.mle <- function(theta.mle,theta){
  bias.mle.mu <- mean(theta.mle[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle <- bias.mle(theta.mle,theta)

save(theta.mle,file="theta.mle.52100.Rdata")

# MOM
mom <- function(x){
  mom.1 <- function(x){  
    sample.x <- sample(x, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x), mom.1(x))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom[i,] <- mom(x[i,])
}
theta.mom.mean <- apply(theta.mom,2,mean)
bias.mom <- function(theta.mom,theta){
  bias.mom.mu <- mean(theta.mom[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom <- bias.mom(theta.mom,theta)

save(theta.mom,file="theta.mom.52100.Rdata")
save(theta.mom.mean,file="theta.mom.mean.52100.Rdata")
save(bias.mom,file="bias.mom.52100.Rdata")

# Bayes
bayes <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes <- bayes(5,2)
theta.bayes.mean <- apply(theta.bayes,2,mean)

bias.bayes <- function(theta.bayes,theta){
  bias.bayes.mu <- mean(theta.bayes[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes <- bias.bayes(theta.bayes,theta)

save(theta.bayes,file="bias.bayes.52100.Rdata")
save(theta.bayes.mean,file="theta.bayes.mean.52100.Rdata")
save(bias.bayes,file="bias.bayes.52100.Rdata")

######################################################################

mu <- 5
beta <- 2
n.1 <- 500
B <- 1000

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x.11 <- replicate(B, Finv(runif(n.1),mu,beta))
x.11 <- matrix(x.11, nrow=B, ncol=n.1, byrow=TRUE)

# MLE
mle.11 <- function(x.11){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon & count<100) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle.11 <- mle.11(x[1,])
theta.mle.11 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle.11[i,] <- mle.11(x[i,])
}
theta.mle.mean.11 <- apply(theta.mle.11,2,mean)
bias.mle.11 <- function(theta.mle.11,theta){
  bias.mle.mu <- mean(theta.mle.11[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle.11[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle.11[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle.11[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle.11 <- bias.mle.11(theta.mle.11,theta)

save(theta.mle.11,file="theta.mle.52500.Rdata")
save(theta.mle.mean.11,file="theta.mle.mean.52500.Rdata")
save(bias.mle.11,file="bias.mle.52500.Rdata")

# MOM
mom.11 <- function(x.11){
  mom.1 <- function(x.11){  
    sample.x <- sample(x.11, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x.11), mom.1(x.11))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom.11 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom.11[i,] <- mom.11(x[i,])
}
theta.mom.mean.11 <- apply(theta.mom.11,2,mean)
bias.mom.11 <- function(theta.mom.11,theta){
  bias.mom.mu <- mean(theta.mom.11[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom.11[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom.11 <- bias.mom.11(theta.mom.11,theta)

save(theta.mom.11,file="theta.mom.52500.Rdata")
save(theta.mom.mean.11,file="theta.mom.mean.52500.Rdata")
save(bias.mom.11,file="bias.mom.52500.Rdata")

# Bayes
bayes.11 <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes.11 <- bayes.11(mu,beta)
theta.bayes.mean.11 <- apply(theta.bayes.11,2,mean)

bias.bayes.11 <- function(theta.bayes.11,theta){
  bias.bayes.mu <- mean(theta.bayes.11[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes.11[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes.11[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes.11[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes.11 <- bias.bayes.11(theta.bayes.11,theta)

save(theta.bayes.11,file="bias.bayes.52500.Rdata")
save(theta.bayes.mean.11,file="theta.bayes.mean.52500.Rdata")
save(bias.bayes.11,file="bias.bayes.52500.Rdata")

######################################################################

mu <- 5
beta <- 2
n.2 <- 1000
B <- 1000

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x.12 <- replicate(B, Finv(runif(n.2),mu,beta))
x.12 <- matrix(x.12, nrow=B, ncol=n.2, byrow=TRUE)

# MLE
mle.12 <- function(x.12){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon & count<100) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle.12 <- mle.12(x[1,])
theta.mle.12 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle.12[i,] <- mle.12(x[i,])
}
theta.mle.mean.12 <- apply(theta.mle.12,2,mean)
bias.mle.12 <- function(theta.mle.12,theta){
  bias.mle.mu <- mean(theta.mle.12[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle.12[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle.12[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle.12[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle.12 <- bias.mle.12(theta.mle.12,theta)

# MOM
mom.12 <- function(x.12){
  mom.1 <- function(x.12){  
    sample.x <- sample(x.12, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x.12), mom.1(x.12))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom.12 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom.12[i,] <- mom.12(x[i,])
}
theta.mom.mean.12 <- apply(theta.mom.12,2,mean)
bias.mom.12 <- function(theta.mom.12,theta){
  bias.mom.mu <- mean(theta.mom.12[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom.12[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom.12 <- bias.mom.12(theta.mom.12,theta)

save(theta.mom.12,file="theta.mom.521000.Rdata")
save(theta.mom.mean.12,file="theta.mom.mean.521000.Rdata")
save(bias.mom.12,file="bias.mom.521000.Rdata")

# Bayes
bayes.12 <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes.12 <- bayes.12(mu,beta)
theta.bayes.mean.12 <- apply(theta.bayes.12,2,mean)

bias.bayes.12 <- function(theta.bayes.12,theta){
  bias.bayes.mu <- mean(theta.bayes.12[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes.12[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes.12[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes.12[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes.12 <- bias.bayes.12(theta.bayes.12,theta)

save(theta.bayes.11,file="bias.bayes.521000.Rdata")
save(theta.bayes.mean.11,file="theta.bayes.mean.521000.Rdata")
save(bias.bayes.11,file="bias.bayes.521000.Rdata")

# Plots
xaxis <- c(100,500,1000)
pdf("sim-parameter-mu-52.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(omi = c(0, 0, .3, 0))
# plot(density(theta.mle[,1]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mom[,1]),col='blue',ylab="Density",xlab="Mean",main="Density of draws of Mean")
par(new=TRUE)
plot(density(theta.bayes[,1]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[1],bias.mle.11[1],bias.mle.12[1]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[1],bias.mom.11[1],bias.mom.12[1]),type="l",col='blue',xaxt = 'n',ylab="Bias",xlab="N",main="Bias of Mean")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[1],bias.bayes.11[1],bias.bayes.12[1]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[3],bias.mle.11[3],bias.mle.12[3]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[3],bias.mom.11[3],bias.mom.12[3]),type="l",col='blue',yaxt='n',ylab="MSE",xlab="N",main="MSE of Mean")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[3],bias.bayes.11[3],bias.bayes.12[3]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
title(main = "Estimator Performance for Location Parameter", outer = TRUE)
dev.off()

pdf("sim-parameter-beta-52.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(omi = c(0, 0, .3, 0))
# plot(density(theta.mle[,2]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mom[,2]),col='blue',ylab="Density",xlab="Scale",main="Density of draws of Scale")
par(new=TRUE)
plot(density(theta.bayes[,2]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[2],bias.mle.11[2],bias.mle.12[2]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[2],bias.mom.11[2],bias.mom.12[2]),type="l",col='blue',xaxt = 'n',ylab="Bias",xlab="N",main="Bias of Scale")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[2],bias.bayes.11[2],bias.bayes.12[2]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[4],bias.mle.11[4],bias.mle.12[4]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[4],bias.mom.11[4],bias.mom.12[4]),type="l",col='blue',ylab="MSe",xlab="N",main="MSE of Scale")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[4],bias.bayes.11[4],bias.bayes.12[4]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
title(main = "Estimator Performance for Location Parameter", outer = TRUE)
dev.off()

######################################################################

mu <- 10
beta <- 4
theta <- c(mu,beta)
n <- 100
B <- 1000

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x <- replicate(B, Finv(runif(n),mu,beta))
x <- matrix(x, nrow=B, ncol=n, byrow=TRUE)

# MLE
mle <- function(x){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle[i,] <- mle(x[i,])
}
theta.mle.mean <- apply(theta.mle,2,mean)
bias.mle <- function(theta.mle,theta){
  bias.mle.mu <- mean(theta.mle[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle <- bias.mle(theta.mle,theta)

save(theta.mle,file="theta.mle.52100.Rdata")

# MOM
mom <- function(x){
  mom.1 <- function(x){  
    sample.x <- sample(x, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x), mom.1(x))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom[i,] <- mom(x[i,])
}
theta.mom.mean <- apply(theta.mom,2,mean)
bias.mom <- function(theta.mom,theta){
  bias.mom.mu <- mean(theta.mom[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom <- bias.mom(theta.mom,theta)

save(theta.mom,file="theta.mom.52100.Rdata")
save(theta.mom.mean,file="theta.mom.mean.52100.Rdata")
save(bias.mom,file="bias.mom.52100.Rdata")

# Bayes
bayes <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes <- bayes(10,4)
theta.bayes.mean <- apply(theta.bayes,2,mean)

bias.bayes <- function(theta.bayes,theta){
  bias.bayes.mu <- mean(theta.bayes[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes <- bias.bayes(theta.bayes,theta)

save(theta.bayes,file="bias.bayes.52100.Rdata")
save(theta.bayes.mean,file="theta.bayes.mean.52100.Rdata")
save(bias.bayes,file="bias.bayes.52100.Rdata")

######################################################################

mu <- 10
beta <- 4
n.1 <- 500
B <- 100

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x.11 <- replicate(B, Finv(runif(n.1),mu,beta))
x.11 <- matrix(x.11, nrow=B, ncol=n.1, byrow=TRUE)

# MLE
mle.11 <- function(x.11){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon & count<100) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}

for(i in 1:B){
  theta.mle.11[i,] <- mle.11(x[i,])
}
theta.mle.mean.11 <- apply(theta.mle.11,2,mean)
bias.mle.11 <- function(theta.mle.11,theta){
  bias.mle.mu <- mean(theta.mle.11[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle.11[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle.11[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle.11[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle.11 <- bias.mle.11(theta.mle.11,theta)

save(theta.mle.11,file="theta.mle.52500.Rdata")
save(theta.mle.mean.11,file="theta.mle.mean.52500.Rdata")
save(bias.mle.11,file="bias.mle.52500.Rdata")

# MOM
mom.11 <- function(x.11){
  mom.1 <- function(x.11){  
    sample.x <- sample(x.11, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x.11), mom.1(x.11))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom.11 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom.11[i,] <- mom.11(x[i,])
}
theta.mom.mean.11 <- apply(theta.mom.11,2,mean)
bias.mom.11 <- function(theta.mom.11,theta){
  bias.mom.mu <- mean(theta.mom.11[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom.11[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom.11 <- bias.mom.11(theta.mom.11,theta)

save(theta.mom.11,file="theta.mom.52500.Rdata")
save(theta.mom.mean.11,file="theta.mom.mean.52500.Rdata")
save(bias.mom.11,file="bias.mom.52500.Rdata")

# Bayes
bayes.11 <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes.11 <- bayes.11(mu,beta)
theta.bayes.mean.11 <- apply(theta.bayes.11,2,mean)

bias.bayes.11 <- function(theta.bayes.11,theta){
  bias.bayes.mu <- mean(theta.bayes.11[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes.11[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes.11[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes.11[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes.11 <- bias.bayes.11(theta.bayes.11,theta)

save(theta.bayes.11,file="bias.bayes.52500.Rdata")
save(theta.bayes.mean.11,file="theta.bayes.mean.52500.Rdata")
save(bias.bayes.11,file="bias.bayes.52500.Rdata")

######################################################################

mu <- 10
beta <- 4
n.2 <- 1000
B <- 1000

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x.12 <- replicate(B, Finv(runif(n.2),mu,beta))
x.12 <- matrix(x.12, nrow=B, ncol=n.2, byrow=TRUE)

# MLE
mle.12 <- function(x.12){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon & count<100) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle.12 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle.12[i,] <- mle.12(x[i,])
}
theta.mle.mean.12 <- apply(theta.mle.12,2,mean)
bias.mle.12 <- function(theta.mle.12,theta){
  bias.mle.mu <- mean(theta.mle.12[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle.12[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle.12[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle.12[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle.12 <- bias.mle.12(theta.mle.12,theta)

# MOM
mom.12 <- function(x.12){
  mom.1 <- function(x.12){  
    sample.x <- sample(x.12, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x.12), mom.1(x.12))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom.12 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom.12[i,] <- mom.12(x[i,])
}
theta.mom.mean.12 <- apply(theta.mom.12,2,mean)
bias.mom.12 <- function(theta.mom.12,theta){
  bias.mom.mu <- mean(theta.mom.12[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom.12[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom.12 <- bias.mom.12(theta.mom.12,theta)

save(theta.mom.12,file="theta.mom.521000.Rdata")
save(theta.mom.mean.12,file="theta.mom.mean.521000.Rdata")
save(bias.mom.12,file="bias.mom.521000.Rdata")

# Bayes
bayes.12 <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes.12 <- bayes.12(mu,beta)
theta.bayes.mean.12 <- apply(theta.bayes.12,2,mean)

bias.bayes.12 <- function(theta.bayes.12,theta){
  bias.bayes.mu <- mean(theta.bayes.12[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes.12[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes.12[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes.12[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes.12 <- bias.bayes.12(theta.bayes.12,theta)

save(theta.bayes.11,file="bias.bayes.521000.Rdata")
save(theta.bayes.mean.11,file="theta.bayes.mean.521000.Rdata")
save(bias.bayes.11,file="bias.bayes.521000.Rdata")

# Plots
xaxis <- c(100,500,1000)
pdf("sim-parameter-mu-104.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(omi = c(0, 0, .3, 0))
# plot(density(theta.mle[,1]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mom[,1]),col='blue',ylab="Density",xlab="Mean",main="Density of draws of Mean")
par(new=TRUE)
plot(density(theta.bayes[,1]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[1],bias.mle.11[1],bias.mle.12[1]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[1],bias.mom.11[1],bias.mom.12[1]),type="l",col='blue',xaxt = 'n',ylab="Bias",xlab="N",main="Bias of Mean")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[1],bias.bayes.11[1],bias.bayes.12[1]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[3],bias.mle.11[3],bias.mle.12[3]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[3],bias.mom.11[3],bias.mom.12[3]),type="l",col='blue',yaxt='n',ylab="MSE",xlab="N",main="MSE of Mean")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[3],bias.bayes.11[3],bias.bayes.12[3]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
title(main = "Estimator Performance for Location Parameter", outer = TRUE)
dev.off()

pdf("sim-parameter-beta-104.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(omi=c(0,0,.3,0))
# plot(density(theta.mle[,2]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mom[,2]),col='blue',ylab="Density",xlab="Scale",main="Density of draws of Scale")
par(new=TRUE)
plot(density(theta.bayes[,2]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[2],bias.mle.11[2],bias.mle.12[2]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[2],bias.mom.11[2],bias.mom.12[2]),type="l",col='blue',xaxt = 'n',ylab="Bias",xlab="N",main="Bias of Scale")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[2],bias.bayes.11[2],bias.bayes.12[2]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[4],bias.mle.11[4],bias.mle.12[4]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[4],bias.mom.11[4],bias.mom.12[4]),type="l",col='blue',ylab="MSe",xlab="N",main="MSE of Scale")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[4],bias.bayes.11[4],bias.bayes.12[4]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
title(main = "Estimator Performance for Scale Parameter", outer = TRUE)
dev.off()


######################################################################

mu <- 100
beta <- 10
theta <- c(mu,beta)
n <- 100
B <- 1000

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x <- replicate(B, Finv(runif(n),mu,beta))
x <- matrix(x, nrow=B, ncol=n, byrow=TRUE)

# MLE
mle <- function(x){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle[i,] <- mle(x[i,])
}
theta.mle.mean <- apply(theta.mle,2,mean)
bias.mle <- function(theta.mle,theta){
  bias.mle.mu <- mean(theta.mle[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle <- bias.mle(theta.mle,theta)

save(theta.mle,file="theta.mle.52100.Rdata")

# MOM
mom <- function(x){
  mom.1 <- function(x){  
    sample.x <- sample(x, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x), mom.1(x))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom[i,] <- mom(x[i,])
}
theta.mom.mean <- apply(theta.mom,2,mean)
bias.mom <- function(theta.mom,theta){
  bias.mom.mu <- mean(theta.mom[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom <- bias.mom(theta.mom,theta)

save(theta.mom,file="theta.mom.52100.Rdata")
save(theta.mom.mean,file="theta.mom.mean.52100.Rdata")
save(bias.mom,file="bias.mom.52100.Rdata")

# Bayes
bayes <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes <- bayes(100,10)
theta.bayes.mean <- apply(theta.bayes,2,mean)

bias.bayes <- function(theta.bayes,theta){
  bias.bayes.mu <- mean(theta.bayes[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes <- bias.bayes(theta.bayes,theta)

save(theta.bayes,file="bias.bayes.52100.Rdata")
save(theta.bayes.mean,file="theta.bayes.mean.52100.Rdata")
save(bias.bayes,file="bias.bayes.52100.Rdata")

######################################################################

mu <- 100
beta <- 10
n.1 <- 500
B <- 100

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x.11 <- replicate(B, Finv(runif(n.1),mu,beta))
x.11 <- matrix(x.11, nrow=B, ncol=n.1, byrow=TRUE)

# MLE
mle.11 <- function(x.11){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon & count<100) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}

for(i in 1:B){
  theta.mle.11[i,] <- mle.11(x[i,])
}
theta.mle.mean.11 <- apply(theta.mle.11,2,mean)
bias.mle.11 <- function(theta.mle.11,theta){
  bias.mle.mu <- mean(theta.mle.11[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle.11[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle.11[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle.11[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle.11 <- bias.mle.11(theta.mle.11,theta)

save(theta.mle.11,file="theta.mle.52500.Rdata")
save(theta.mle.mean.11,file="theta.mle.mean.52500.Rdata")
save(bias.mle.11,file="bias.mle.52500.Rdata")

# MOM
mom.11 <- function(x.11){
  mom.1 <- function(x.11){  
    sample.x <- sample(x.11, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x.11), mom.1(x.11))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom.11 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom.11[i,] <- mom.11(x[i,])
}
theta.mom.mean.11 <- apply(theta.mom.11,2,mean)
bias.mom.11 <- function(theta.mom.11,theta){
  bias.mom.mu <- mean(theta.mom.11[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom.11[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom.11 <- bias.mom.11(theta.mom.11,theta)

save(theta.mom.11,file="theta.mom.52500.Rdata")
save(theta.mom.mean.11,file="theta.mom.mean.52500.Rdata")
save(bias.mom.11,file="bias.mom.52500.Rdata")

# Bayes
bayes.11 <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes.11 <- bayes.11(mu,beta)
theta.bayes.mean.11 <- apply(theta.bayes.11,2,mean)

bias.bayes.11 <- function(theta.bayes.11,theta){
  bias.bayes.mu <- mean(theta.bayes.11[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes.11[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes.11[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes.11[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes.11 <- bias.bayes.11(theta.bayes.11,theta)

save(theta.bayes.11,file="bias.bayes.52500.Rdata")
save(theta.bayes.mean.11,file="theta.bayes.mean.52500.Rdata")
save(bias.bayes.11,file="bias.bayes.52500.Rdata")

######################################################################


mu <- 100
beta <- 10
n.2 <- 1000
B <- 1000

# Simulation
Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
x.12 <- replicate(B, Finv(runif(n.2),mu,beta))
x.12 <- matrix(x.12, nrow=B, ncol=n.2, byrow=TRUE)

# MLE
mle.12 <- function(x.12){
  make.functions <- function(x) {
    
    n <- length(x)
    
    # Log of the likelihood
    log.likelihood <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z)))
    }
    
    # First derivative of the log likelihood
    gradient <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      c(n/beta - 2*sum((exp(z)/(beta*(exp(z)+1)))), 
        sum((x-mu)/beta^2) - n/beta - 2*sum((x-mu)*exp(z)/(beta^2*(exp(z)+1))))
    }
    
    # Second derivative of the log likelihood
    hessian <- function(theta) {
      mu <- theta[1]
      beta <- theta[2]
      z <- -(x-mu)/beta
      matrix(
        c(-2*sum(exp(z)/(beta^2*(1+exp(z))^2)), 
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -n/beta^2 + 2*sum(exp(z)/(beta^2*(exp(z)+1))) - 2*sum((x-mu)*exp(z)/(beta^3*(exp(z)+1))) + 2*sum((x-mu)*exp(2*z)/(beta^3*(exp(z)+1)^2)),
          -2*sum((x-mu)/beta^3) + n/beta^2 - 2*sum((x-mu)^2*exp(z)/(beta^2*(1+exp(z)))) + 2*sum((x-mu)^2*exp(2*z)/(beta^2*(exp(z)+1)^2))),
        nrow=2)
    }
    
    list(logLike=log.likelihood,gradient=gradient,hessian=hessian)
    
  }  
  # function to randomly pull from a logistic model
  f <- make.functions(x)
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon & count<100) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle.12 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle.12[i,] <- mle.12(x[i,])
}
theta.mle.mean.12 <- apply(theta.mle.12,2,mean)
bias.mle.12 <- function(theta.mle.12,theta){
  bias.mle.mu <- mean(theta.mle.12[,1] - theta[1])
  bias.mle.beta <- mean(theta.mle.12[,2] - theta[2])
  mse.mle.mu <- bias.mle.mu^2 + var(theta.mle.12[,1])
  mse.mle.beta <- bias.mle.beta^2 + var(theta.mle.12[,2])
  return(c(bias.mle.mu, bias.mle.beta, mse.mle.mu, mse.mle.beta))
}
bias.mle.12 <- bias.mle.12(theta.mle.12,theta)

# MOM
mom.12 <- function(x.12){
  mom.1 <- function(x.12){  
    sample.x <- sample(x.12, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x.12), mom.1(x.12))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}
theta.mom.12 <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom.12[i,] <- mom.12(x[i,])
}
theta.mom.mean.12 <- apply(theta.mom.12,2,mean)
bias.mom.12 <- function(theta.mom.12,theta){
  bias.mom.mu <- mean(theta.mom.12[,1] - theta[1])
  bias.mom.beta <- mean(theta.mom.12[,2] - theta[2])
  mse.mom.mu <- bias.mom.mu^2 + var(theta.mom[,1])
  mse.mom.beta <- bias.mom.beta^2 + var(theta.mom[,2])
  return(c(bias.mom.mu, bias.mom.beta, mse.mom.mu, mse.mom.beta))
}
bias.mom.12 <- bias.mom.12(theta.mom.12,theta)

save(theta.mom.12,file="theta.mom.521000.Rdata")
save(theta.mom.mean.12,file="theta.mom.mean.521000.Rdata")
save(bias.mom.12,file="bias.mom.521000.Rdata")

# Bayes
bayes.12 <- function(mu,beta){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) - mu
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) - beta
  }
  Finv <- function(y,mu,beta) -log(1/y-1)*beta + mu
  x <- Finv(runif(n),mu,beta)
  # set initial values
  n.draws <- 10000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  draws[1,1] <- 10
  draws[1,2] <- 10
  prop.sd.mu <- 1
  prop.sd.beta <- 1
  accept.mu <- 0
  accept.beta <- 0
  # Gibbs sampling for both the mu and beta
  for(i in 2:nrow(draws)){
    mu.current <- draws[i-1,1]
    beta.current <- draws[i-1,2]
    mu.proposal <- rnorm(1,mu.current,prop.sd.mu)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    beta.proposal <- rnorm(1,beta.current,prop.sd.beta)
    methast.mu <- complete.mu(x, mu.proposal,beta.current) - complete.mu(x,mu.current,beta.current)
    if(log(runif(1)) < methast.mu){
      draws[i,1] <- mu.proposal
      accept.mu <- accept.mu + 1
    } else{
      draws[i,1] <- mu.current
    }
    if ( beta.proposal > 0 ) {
      methast.beta <- complete.beta(x,mu.current,beta.proposal) - complete.beta(x,mu.current,beta.current)
      if(log(runif(1))<methast.beta){
        draws[i,2] <- beta.proposal
        accept.beta <- accept.beta + 1
      } else{
        draws[i,2] <- beta.current
      }
    } else {
      draws[i,2] <- beta.current
    }
  } 
  return(draws)
}
theta.bayes.12 <- bayes.12(mu,beta)
theta.bayes.mean.12 <- apply(theta.bayes.12,2,mean)

bias.bayes.12 <- function(theta.bayes.12,theta){
  bias.bayes.mu <- mean(theta.bayes.12[,1] - theta[1])
  bias.bayes.beta <- mean(theta.bayes.12[,2] - theta[2])
  mse.bayes.mu <- bias.bayes.mu^2 + var(theta.bayes.12[,1])
  mse.bayes.beta <- bias.bayes.beta^2 + var(theta.bayes.12[,2])
  return(c(bias.bayes.mu, bias.bayes.beta, mse.bayes.mu, mse.bayes.beta))
}
bias.bayes.12 <- bias.bayes.12(theta.bayes.12,theta)

save(theta.bayes.11,file="bias.bayes.521000.Rdata")
save(theta.bayes.mean.11,file="theta.bayes.mean.521000.Rdata")
save(bias.bayes.11,file="bias.bayes.521000.Rdata")

# Plots
xaxis <- c(100,500,1000)
pdf("sim-parameter-mu-10010.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(omi=c(0,0,.3,0))
# plot(density(theta.mle[,1]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mom[,1]),col='blue',ylab="Density",xlab="Mean",main="Density of draws of Mean")
par(new=TRUE)
plot(density(theta.bayes[,1]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[1],bias.mle.11[1],bias.mle.12[1]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[1],bias.mom.11[1],bias.mom.12[1]),type="l",col='blue',xaxt = 'n',ylab="Bias",xlab="N",main="Bias of Mean")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[1],bias.bayes.11[1],bias.bayes.12[1]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[3],bias.mle.11[3],bias.mle.12[3]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[3],bias.mom.11[3],bias.mom.12[3]),type="l",col='blue',yaxt='n',ylab="MSE",xlab="N",main="MSE of Mean")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[3],bias.bayes.11[3],bias.bayes.12[3]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
title(main = "Estimator Performance for Location Parameter", outer = TRUE)
dev.off()

pdf("sim-parameter-beta-10010.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(omi=c(0,0,.3,0))
# plot(density(theta.mle[,2]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mom[,2]),col='blue',ylab="Density",xlab="Scale",main="Density of draws of Scale")
par(new=TRUE)
plot(density(theta.bayes[,2]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[2],bias.mle.11[2],bias.mle.12[2]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[2],bias.mom.11[2],bias.mom.12[2]),type="l",col='blue',xaxt = 'n',ylab="Bias",xlab="N",main="Bias of Scale")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[2],bias.bayes.11[2],bias.bayes.12[2]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')

# plot(xaxis,c(bias.mle[4],bias.mle.11[4],bias.mle.12[4]),col='purple',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
# par(new=TRUE)
plot(xaxis,c(bias.mom[4],bias.mom.11[4],bias.mom.12[4]),type="l",col='blue',ylab="MSe",xlab="N",main="MSE of Scale")
ticks<-c(100,500,1000)
axis(1,at=ticks,labels=ticks)
par(new=TRUE)
plot(xaxis,c(bias.bayes[4],bias.bayes.11[4],bias.bayes.12[4]),type="l",col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
title(main = "Estimator Performance for Scale Parameter", outer = TRUE)
dev.off()