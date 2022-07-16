####################################################################
########## Parameter Estimation for Beta Regression Model ##########
####################################################################

# Estimation using the Quasi-Newton method with optim()

EnvWomen <- read.xls(("environment_data.xlsx"),sheet=1,header=TRUE)

## install.packages("betareg")
library(betareg)

data("EnvWomen", package="betareg")
data_logit <- betareg(y~x1+x2+x3+x4, data=dat,link="logit", method="BFGS")
summary(data_logit)


# Estimation using the Newton-Raphson method: Decided not to use, but these are the general steps for how to perform the MLE approximation using this method

##################################
######### Link function ##########
##################################

# logit link function
g <- function(y) sapply(y,function(y) log(y/(1-y)))

# g'(mu_t)
dg <- function(mu) sapply(mu, function(mu) -1/((mu-1)*mu))

##################################
##### Log-Likelihood function ####
##################################

ll <- function(mu,phi,n,y) {
  lt <- rep(NA,n)
  for(i in 1:n) {
  lt[i] <- log(gamma(phi)) - log(gamma(mu[i]*phi)) - log(gamma((1-mu[i])*phi)) +
        (mu[i] - 1)*log(y[i]) + ((1 - mu[i])*phi - 1)*log(1-y[i])
  }
  ll <- sum(lt)
  return(ll)
}

#################################
####### Gradient function #######
#################################

gradient.ll <- function(mu,phi,n,y,X,dg=dg) {
  y.star <- sapply(y,function(y) log(y/(1-y)))
  #mu.star <- sapply(mu, function(mu) digamma(mu*phi)-digamma((1-mu)*phi))
  mu.star <- rep(NA,n)
  for(i in 1:n){
    mu.star[i] <- digamma(mu[i]*phi)-digamma((1-mu[i])*phi)
  }
  
  x <- rep(NA,n)
  dg <- dg(mu)
  for(i in 1:n) x[i] <- 1/dg[i]
  T.diag <- diag(x=x)
  
  # Score function for beta
  U.beta <- phi*t(X)%*%T.diag%*%(y.star - mu.star)
  
  # Score function for phi
  #U.phi.i <- rep(NA,n)
  #for(i in 1:n)
  #  U.phi.i[i] <- mu[i]*(y.star[i] - mu.star[i]) + log(1-y[i]) - 
  #    digamma((1-mu[i])*phi) + digamma(phi)
  #U.phi <- sum(U.phi.i)
  
  return(matrix(c(t(U.beta)),ncol=1))
}

###################################
##### Hessian Matrix Function #####
###################################

# Set up functions for the hessian matrix
w <- function(y,X,mu,phi,n,dg) {
  w <- rep(NA,n)
  dgmu <- dg(mu)
  for(i in 1:n) {
    w[i] <- phi*(trigamma(mu[i]*phi) + trigamma((1-mu[i])*phi))/(dgmu[i]^2)
  }
  return(diag(x=w))
}

cf <- function(mu,phi,n) {
  c <- rep(NA,n)
  for(i in 1:n) {
    c[i] <- phi*(trigamma(mu[i]*phi)*mu[i] - trigamma((1-mu[i])*phi)*(1-mu[i]))
  }
  return(t(c))
}

d <- function(mu,phi,n) {
  d <- rep(NA,n)
  for(i in 1:n) {
    d[i] <- trigamma(mu[i]*phi)*(mu[i])^2 + trigamma((1-mu[i])*phi)*(1-mu[i])^2 - trigamma(phi)
  }
  return(diag(x=d))  
}

# Hessian matrix:
hessian.ll <- function(mu,phi,X,w,cf,d,dg) {
  W <- w(y,X,mu,phi,n,dg)
  c <- cf(mu,phi,n)
  D <- d(mu,phi,n)
  
  x <- rep(NA,n)
  dgmu <- dg(mu)
  for(i in 1:n) x[i] <- 1/dgmu[i]
  T.diag <- diag(x=x)
  
  K.bb <- phi*t(X)%*%W%*%X
  K.bp <- t(X)%*%T.diag%*%t(c)
  K.pb <- t(K.bp)
  K.pp <- sum(D)
  return(matrix(c(K.bb,K.bp,K.pb,K.pp),byrow=T,nrow=6,ncol=6))
}

# Inverse of hessian matrix: variance-covariance matrix for multivariate normal distribution of beta and phi
info.matrix.inv <- function(mu,phi,n,X,W,c) {
  T.diag <- diag(x=sapply(mu, 1/dg(mu)))
  I <- diag(n)
  gam <- tr(D) - solve(phi)%*%t(c)%*%t(T.diag)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%T.diag%*%c
  K.bb <- 1/phi*solve(t(X)%*%W%*%X)%*%(I+(t(X)%*%T.diag%*%c%*%t(c)%*%t(T.diag)%*%X%*%solve(t(X)%*%W%*%X))/(gam*phi))
  K.bp <- 1/(gam*phi)%*%solve(t(X)%*%W%*%X)%*%t(X)%*%T.diag%*%c
  K.pb <- t(K.bp)
  K.pp <- solve(gam)
  return(matrix(c(K.bb,K.bp,K.pb,K.pp),byrow=T,nrow=2,ncol=2))
}

############################################
### Simulation to find mle of parameters ###
############################################

# find initial values
beta.init <- function(X,y,ftn=g){
  z <- ftn(y)
  return(solve(t(X)%*%X)%*%t(X)%*%z)
}

phi.init <- function(X,y,n,k,ftn1=g,ftn2=dg){
  z <- ftn1(y)
  xbeta <- rep(NA,n)
  for(i in 1:n) {
    xbeta[i] <- X[i,]%*%solve(t(X)%*%X)%*%t(X)%*%z
  }
  mu <- sapply(xbeta, function(y=xbeta) exp(y)/(1+exp(y)))
  
  dgmu <- ftn2(mu)
  e <- z - X%*%solve(t(X)%*%X)%*%t(X)%*%z
  
  var.y <- rep(NA,n)
  phi.t <- rep(NA,n)
  for(t in 1:n) {
    var.y[t] <- t(e)%*%e/((n-k)*(dgmu[t])^2)
    phi.t[t] <- mu[t]*(1-mu[t])/var.y[t]
  }
  return(1/n*sum(phi.t) - 1)
}
 
phi <- .5 
theta <- matrix(beta.init(X,y),ncol=1)

mu.init <- function(beta,X,n) {
  mu.init <- rep(NA,n)
  for(i in 1:n) mu.init[i] <- exp(t(X[i,])%*%beta)/(1+exp(t(X[i,])%*%beta))
  return(mu.init)
}
  
epsilon <- 0.0001
count <- 0
# find the mle
niter <- 10
while(count < niter){
while ( sqrt(sum(gradient.ll(mu.init(theta,X,n),phi,n,y,X,dg)^2)) > epsilon) {
  theta <- theta - solve(hessian.ll(mu.init(theta,X,n),phi,X,w,cf,d,dg)[-6,-6], matrix(gradient.ll(mu.init(theta,X,n),phi,n,y,X,dg),nrow=k))
  count <- count + 1
}
}

