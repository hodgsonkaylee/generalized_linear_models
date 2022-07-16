###################
### Import Data ###
###################

# install.packages("gdata")
library(gdata)

EnvWomen <- read.xls(("environment_data.xlsx"),sheet=1,header=TRUE)

#################################
### Estimating the Parameters ###
#################################

## install.packages("betareg")
library(betareg)

# Beta regression fit
data_logit <- betareg(EnvDeaths~womeninparliament+genderequalityfw+govtparticipationscale+womenministerial, data=EnvWomen, 
                      link="logit", method="BFGS")
# Parameter Estimates
paramest <- coef(data_logit)
# Mean Squared Error
mean(data_logit$residuals^2)

###################################
### Hypothesis testing for Data ###
###################################

# Reduced Models:
data_logit1 <- betareg(EnvDeaths~genderequalityfw+govtparticipationscale+womenministerial, data=EnvWomen, 
                      link="logit", method="BFGS")
data_logit2 <- betareg(EnvDeaths~womeninparliament+govtparticipationscale+womenministerial, data=EnvWomen, 
                      link="logit", method="BFGS")
data_logit3 <- betareg(EnvDeaths~womeninparliament+genderequalityfw+womenministerial, data=EnvWomen, 
                      link="logit", method="BFGS")
data_logit4 <- betareg(EnvDeaths~womeninparliament+genderequalityfw+govtparticipationscale, data=EnvWomen, 
                      link="logit", method="BFGS")

library(lmtest)

# Wald Test
waldtest(data_logit,data_logit1,test="F")
waldtest(data_logit,data_logit2,test="F")
waldtest(data_logit,data_logit3,test="F")
waldtest(data_logit,data_logit4,test="F")


##############################
### Predicting Probability ###
##############################

# new observations of the variables
Xnewh <- c(1,0.39,0,0,0.39)
Xnewm <- c(1,0.25,3,2,0.25)
Xnewl <- c(1,0.10,5,4,0.10)
phi <- paramest[6]
# calculate mu from the new observations
mu <- exp(t(Xnew)%*%paramest[1:5])/(1+exp(t(Xnew)%*%paramest[1:5]))
# beta pdf function
f <- function(y) gamma(phi)/(gamma(mu*phi)*gamma((1-mu)*phi))*y^(mu*phi-1)*(1-y)^((1-mu)*phi-1)

# Use Simpson's rule to approximate integral of pdf
simpson <- function(ftn=f, a=0, b=k, n=100,mu=mu,phi=phi,less=F){
  n <- max(c(2*(n%/%2),4))
  h <- (b-a)/n
  x.vec1 <- seq(a+h, b-h, by=2*h)
  x.vec2 <- seq(a+2*h, b-2*h, by=2*h)
  f.vec1 <- sapply(x.vec1,ftn)
  f.vec2 <- sapply(x.vec2,ftn)
  S <- h/3*(ftn(a) + ftn(b) + 4*sum(f.vec1) + 2*sum(f.vec2))
  int <- integrate(f,a,b)$val
  error <- abs(S-int)
  if(less==F) return(c(1-S,error))
  else if(less==T) return(c(S,error))
}

# Find the probability that y>0.20
k <- 0.20
mu <- exp(t(Xnewl)%*%paramest[1:5])/(1+exp(t(Xnewl)%*%paramest[1:5]))
problow <- simpson()
mu <- exp(t(Xnewm)%*%paramest[1:5])/(1+exp(t(Xnewm)%*%paramest[1:5]))
probmed <- simpson()
mu <- exp(t(Xnewh)%*%paramest[1:5])/(1+exp(t(Xnewh)%*%paramest[1:5]))
probhigh <- simpson()

# Find the probability that y<0.10
k <- 0.10
mu <- exp(t(Xnewl)%*%paramest[1:5])/(1+exp(t(Xnewl)%*%paramest[1:5]))
problow <- simpson(less=T)
mu <- exp(t(Xnewm)%*%paramest[1:5])/(1+exp(t(Xnewm)%*%paramest[1:5]))
probmed <- simpson(less=T)
mu <- exp(t(Xnewh)%*%paramest[1:5])/(1+exp(t(Xnewh)%*%paramest[1:5]))
probhigh <- simpson(less=T)


##############################
### Linear Fit of the Data ###
##############################

# Find the beta coefficients for a binomial fit of the data with the logit link function
fit.dat.b <- glm(EnvDeaths~womeninparliament+genderequalityfw+govtparticipationscale+womenministerial, data=EnvWomen, family=quasibinomial(link = "logit"))
# parameter estimates
summary(fit.dat)
# Mean Square error
mean(weighted.residuals(fit.dat.b)^2)




