## Contraceptive Data

####################################
  # Import Contraceptive Data #
####################################

data.cont <- c(23.00,69.00,57.00,18.00,77.00,55.00,70.00,55.00,62.00,59.00,63.00,70.00,55.00,18.00,66.00,61.00,46.00,
           53.00,69.00,17.00,46.00,22.00,56.00,34.00,15.00,6.00,88.00,79.00,19.00,45.00,76.00,18.00,74.00,86.00,
           20.00,19.00,70.00,22.00,80.00,59.00,72.00,13.00,8.00,34.00,44.00,76.00,31.00,9.00,53.00,7.00,61.00,6.00,
           16.00,34.00,35.00,73.00,55.00,63.00,77.00,53.00,73.00,61.00,51.00,58.00,42.00,50.00,54.00,60.00,20.00,
           42.00,40.00,40.00,59.00,35.00,10.00,11.00,73.00,60.00,55.00,23.00,67.00,12.00,56.00,50.00,69.00,80.00,
           14.00,15.00,71.00,30.00,35.00,57.00,63.00,32.00,79.00,75.00,55.00,38.00,68.00,53.00,24.00,22.00,58.00,
           17.00,35.00,80.00,4.00,68.00,12.00,48.00,66.00,54.00,28.00,38.00,79.00,20.00,63.00,74.00,27.00,65.00,
           84.00,76.00,49.00,76.00,34.00,49.00,67.00)

####################################
  # Plot Histogram and Density #
####################################

pdf("Data-Density.pdf")
par(mfrow=c(1,1),cex.axis=1,cex.lab=1.5,cex.main=1.5)
x <- seq(0,100,length=100)
hist(data.cont, freq=NULL,ylab="Density",xlab="Access to Contraceptives (%)", main="Contraceptive Data Density Plot")
par(new=TRUE)
plot(density(data.cont),type="l",col="blue",bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(x,dlogis(x,mu,beta),type="l",col="green",bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
legend("topright",c("Data","Logistic"), col=c("blue","green"), lwd=c(1,1))
dev.off()

n <- length(data.cont)

####################################
    # Sample from the Data #
####################################

B <- 1000
data <- replicate(B,sample(data.cont,replace=T))
data <- matrix(data, nrow=B, ncol=length(data.cont), byrow=TRUE)

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
  f <- make.functions(x)
  # initial values
  theta <- c(40,10)
  epsilon <- 0.0001
  count <- 0
  # find the mle
  while ( sqrt(sum(f$gradient(theta)^2)) > epsilon ) {
    theta <- theta - solve(f$hessian(theta), f$gradient(theta))
    count <- count + 1
  }
  return(theta)
}
theta.mle.data <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mle.data[i,] <- mle(data[i,])
}

# Bootstrap confidence interval
theta.mle.data.mean <- apply(theta.mle.data,2,mean)
ci.mle.data <- apply(theta.mle.data,2,quantile,c(.025, .975))

####################################
      # Method of Moments #
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
theta.mom.data <- matrix(0,ncol=2,nrow=B)
for(i in 1:B){
  theta.mom.data[i,] <- mom(data[i,])
}

# bootstrap confidence interval
theta.mom.data.mean <- apply(theta.mom.data,2,mean)
ci.mom.data <- apply(theta.mom.data,2,quantile,c(.025, .975))

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
    -sum((x-mu)/beta) - n*log(beta) - 2*sum(log(1+exp(z))) + log(1/sqrt(2*pi)) - beta/2
  }
  # set initial values
  n.draws <- 12000
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
  mu.hat <- mean(draws[,1])
  beta.hat <- mean(draws[,2])
  return(draws[2001:12000,])
}
theta.bayes.data <- bayes(data.cont)

# credible interval for the draws
theta.bayes.data.mean <- apply(theta.bayes.data, 2, mean)
ci.bayes.data <- apply(theta.bayes.data,2,quantile,c(.025,.975))

########################################
# Distribution of Parameter Estimates #
########################################

pdf("data-parameter-est.pdf",width=12,height=6)
par(mfrow=c(1,2))
par(cex.axis=1)
par(cex.lab=1.5)
par(cex.main=2)
par(mgp = c(1.8, .7, 0))
plot(density(theta.mom.data[,1]),col='blue',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mle.data[,1]),col='purple',ylab="Density",xlab="Mean",main="Density of draws of Mean")
par(new=TRUE)
plot(density(theta.bayes.data[,1]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
legend("topleft",c("MOM","MLE","Bayes"), cex=1.1, col=c("blue","purple","red"), lwd=c(1,1,1))

plot(density(theta.mom.data[,2]),col='blue',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
par(new=TRUE)
plot(density(theta.mle.data[,2]),col='purple',ylab="Density",xlab="Mean",main="Density of draws of Beta")
par(new=TRUE)
plot(density(theta.bayes.data[,2]),col='red',bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
dev.off()

######################################
# Density Curves based on Estimates #
######################################
pdf("data-distribution-est.pdf")
par(mfrow=c(1,1),cex.axis=1,cex.lab=1.5,cex.main=1.5)
x <- seq(0,100,length=100)
y <- seq(0,0.02,length=100)
plot(x,y,type="n",ylab="Density",xlab="Contraceptive Data",main="Density Curves from the Estimates",ylim = c(0,0.02))
points(density(data.cont,from=0,to=100),type="l",col="black",bty = 'n', ann = FALSE, xaxt = 'n', yaxt = 'n')
lines(x, dlogis(x,theta.mom.data.mean[1],theta.mom.data.mean[2]), lty=2, col = "blue")
lines(x, dlogis(x,theta.mle.data.mean[1],theta.mle.data.mean[2]),lty=2,col = "purple")
lines(x, dlogis(x,theta.bayes.data.mean[1],theta.bayes.data.mean[2]), lty=2, col = "red")
legend("topleft",c("Dataset","MOM","MLE","Bayes"), 
       cex=1.2, col=c("black","blue","purple","red"), lwd=c(1,1,1,1))
dev.off()
