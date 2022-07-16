#####################################
# Model Mispecification Simulation #
#####################################

# Random normal generator through box-muller method
rnorm <- function(X, M=0, V=1){
  rnorms <- function(M=0, V=1){    
    u1 <- runif(1)
    u2 <- runif(1)
    x1 <- M + V*sqrt(-2*log(u1))*cos(2*pi*u2)
    x2 <- M + V*sqrt(-2*log(u1))*sin(2*pi*u2)
    return(x1)
  }
  x <- replicate(X, rnorms(M,V))
  return(x)
}

####################################
# Maximum Likelihood Estimators #
####################################

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
  # initial values
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  f <- make.functions(x)
  if (theta[2]<1){
    stop("Starting value must be greater than zero.")
  }
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon) {
    theta <- theta - solve(f$hessian(theta), matrix(f$gradient(theta)))
    count <- count + 1
  }
  return(theta)
}

####################################
# Method of Moments Estimators #
####################################

mom <- function(x){
  mom.1 <- function(x){  
    sample.x <- sample(x, replace=TRUE) # sample from the sample data instead of the population
  }
  mom.rep <- replicate(length(x), mom.1(x))
  mu.hat <- mean(mom.rep)
  beta.hat <- sqrt(3*sd(mom.rep)^2/pi^2)
  return(c(mu.hat,beta.hat))
}

####################################
# Bayesian Estimators #
####################################

bayes <- function(x){
  complete.mu <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - 2*sum(log(1+exp(z))) + log(1/sqrt(2*pi)) - mu/2
  }
  complete.beta <- function(x,mu,beta){
    n <- length(x)
    z <- -(x-mu)/beta
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) + log(1/sqrt(2*pi)) - mu/2
  }
  # set initial values
  n.draws <- 12000
  draws <- matrix(NA,nrow=n.draws,ncol=2)
  mom <- function(x){
    mu.hat <- mean(x)
    beta.hat <- sqrt(3*sd(x)^2/pi^2)
    return(c(mu.hat,beta.hat))
  }
  theta <- mom(x)
  draws[1,1] <- theta[1]
  draws[1,2] <- theta[2]
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
  return(draws[2001:12000,])
}

####################################
      # Bias Estimation #
####################################

mean.estimation.mis <- function(n.runs = 1000, funct, mu, sigma, n = c(100,500,1000)){
  n.Bs <- 2*length(n)
  mean.vec <- matrix(0,n.Bs,n.runs)
  
  for(m in 0:(n.Bs/2-1)){
    for(l in 1:n.runs){
      p.vals <- funct(rnorm(n[m+1], mu, sigma))
      while(abs(p.vals[1]) > 1e4 | abs(p.vals[2]) > 1e10){
        p.vals <- funct(rnorm(n[m+1], mu, sigma))
      }
      mean.vec[(2*(m)+1),l] <- p.vals[1]
      mean.vec[(2*(m+1)),l] <- p.vals[2]
    }
    print(m)
  }
  mean.vec
}
mean.estimation.bayes.mis <- function(funct,mu,sigma,n=c(100,500,1000)){
  n.Bs <- 2*length(n)
  mean.vec <- matrix(0,nrow=n.Bs,ncol=10000)
  for(m in 0:(n.Bs/2-1)){
    p.vals <- funct(rnorm(n[m+1], mu, sigma))
    while(abs(p.vals[1]) > 1e4 | abs(p.vals[2]) > 1e10){
      p.vals <- funct(rnorm(n[m+1], mu, sigma))
    }
    mean.vec[(2*m+1),] <- p.vals[,1]
    mean.vec[(2*(m+1)),] <- p.vals[,2]
    print(m)
  }
  mean.vec
}

mean.mom.mm <- mean.estimation.mis(funct=mom,mu=50,sigma=25)
mean.mle.mm <- mean.estimation.mis(funct=mle,mu=50,sigma=25)
mean.bayes.mm <- mean.estimation.bayes.mis(funct=bayes,mu=50,sigma=25)

############################################
        # Bias Estimation Mean #
############################################

mean.estimates.mis <- function(mean, n = c(100,500,1000)){
  n.val <- dim(mean)[1]/2
  mu.mean <- apply(mean[2*(0:(n.val-1))+1,],1,mean)
  beta.mean <- apply(mean[2*(1:n.val),],1,mean)
  mean.est <- data.frame(n, mu.mean, beta.mean)
  mean.est
}

mom.mean.mm <- mean.estimates.mis(mean.mom.mm)
mle.mean.mm <- mean.estimates.mis(mean.mle.mm)
bayes.mean.mm <- mean.estimates.mis(mean.bayes.mm)

######################################
# Density Curves based on Estimates #
######################################

pdf("data-distribution-mm.pdf")
par(mfrow=c(1,1),cex.axis=1,cex.lab=1.5,cex.main=1.5)
x <- seq(0,100,length=100)
y <- seq(0,0.02,length=100)
plot(x,y,type="n",ylab="Model Mispecification Data",xlab="Density",main="Density Curves of the Model Mispecification Estimates",
     ylim = c(0,0.02))
points(density(data.cont,from=0,to=100),type="l",col="black",bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
lines(x, dlogis(x,mom.mean.mm[1,2],mom.mean.mm[1,3]), lty=1, col = "blue")
lines(x, dlogis(x,mle.mean.mm[1,2],mle.mean.mm[1,3]),lty=1,col = "purple")
lines(x, dlogis(x,bayes.mean.mm[1,2],bayes.mean.mm[1,3]), lty=1, col = "red")
lines(x, dlogis(x,mom.mean.mm[2,2],mom.mean.mm[2,3]), lty=3, col = "blue")
lines(x, dlogis(x,mle.mean.mm[2,2],mle.mean.mm[2,3]),lty=3,col = "purple")
lines(x, dlogis(x,bayes.mean.mm[2,2],bayes.mean.mm[2,3]), lty=3, col = "red")
lines(x, dlogis(x,mom.mean.mm[3,2],mom.mean.mm[3,3]), lty=4, col = "blue")
lines(x, dlogis(x,mle.mean.mm[3,2],mle.mean.mm[3,3]),lty=4,col = "purple")
lines(x, dlogis(x,bayes.mean.mm[3,2],bayes.mean.mm[3,3]), lty=4, col = "red")
legend("topleft",c("Data","MOM","MLE","Bayes","n=100","n=500","n=1000"), 
       cex=1.2, col=c("black","blue","purple","red","black","black","black"), lwd=c(1,1,1,1,1,1,1),lty=c(1,1,1,1,1,3,4))
dev.off()
