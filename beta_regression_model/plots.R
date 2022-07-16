library(ggplot2)
library(grid)
library(gridExtra)

######################################
######### Beta Density Plot ##########
######################################

nsim <- 900000
# p=.5,q=.5, p=2,q=2, p=1,q=5, p=5,q=1
dens.data <- data.frame(parameters = factor(rep(c("p=q=0.5","p=q=2","p=1, q=5","p=5, q=1"), each=nsim)),y=c(rbeta(nsim,.5,.5),rbeta(nsim,2,2),rbeta(nsim,1,5),rbeta(nsim,5,1)))

pdf("betadensity.pdf")

ggplot(dens.data, aes(x=y,colour=parameters)) + geom_density() +
  theme(legend.position = c(0.7, 0.85))

dev.off()

#############################################
######### Simulated y Density Plot ##########
#############################################

simdat1 <- data.frame(Sample = factor(rep(c("n=150","n=1200","n=6000"), c(150,1200,6000))),y=c(sim11$X6,sim12$X6,sim13$X6))
simdat2 <- data.frame(Sample = factor(rep(c("n=150","n=1200","n=6000"), c(150,1200,6000))),y=c(sim21$X6,sim22$X6,sim23$X6))
simdat3 <- data.frame(Sample = factor(rep(c("n=150","n=1200","n=6000"), c(150,1200,6000))),y=c(sim31$X6,sim32$X6,sim33$X6))
simdat4 <- data.frame(Sample = factor(rep(c("n=150","n=1200","n=6000"), c(150,1200,6000))),y=c(sim41$X6,sim42$X6,sim43$X6))


pdf("ydensity.pdf")

simplot1 <- ggplot(simdat1, aes(x=y,colour=Sample)) + geom_density() + labs(x="simulated y") +
  ggtitle(expression(paste("Density of y for ", phi, "=2"))) + theme(legend.position="none") + xlim(0,1)
simplot2 <- ggplot(simdat2, aes(x=y,colour=Sample)) + geom_density() + labs(x="simulated y") +
  ggtitle(expression(paste("Density of y for ", phi, "=5"))) + theme(legend.position="none") + xlim(0,1)
simplot3 <- ggplot(simdat3, aes(x=y,colour=Sample)) + geom_density() + labs(x="simulated y") + 
  ggtitle(expression(paste("Density of y for ", phi, "=10"))) + theme(legend.position="none") + xlim(0,1)
simplot4 <- ggplot(simdat4, aes(x=y,colour=Sample)) + geom_density() + labs(x="simulated y") +
  ggtitle(expression(paste("Density of y for ", phi, "=60"))) + theme(legend.position = c(0.8, 0.7)) + xlim(0,1)
grid.arrange(simplot1,simplot2,simplot3,simplot4,ncol=2)

dev.off()

#########################################
######### Absolute Error Plots ##########
#########################################

erdat1 <- data.frame(Coefficient = factor(c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5)),
                     y = c(sim11$X7[1:6],sim12$X7[1:6],sim13$X7[1:6]),
                     x = factor(rep(c(150,1200,6000),each=6)))
erdat2 <- data.frame(Coefficient = factor(c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5)),
                     y = c(sim21$X7[1:6],sim22$X7[1:6],sim23$X7[1:6]),
                     x = factor(rep(c(150,1200,6000),each=6)))
erdat3 <- data.frame(Coefficient = factor(c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5)),
                     y = c(sim31$X7[1:6],sim32$X7[1:6],sim33$X7[1:6]),
                     x = factor(rep(c(150,1200,6000),each=6)))
erdat4 <- data.frame(Coefficient = factor(c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5)),
                     y = c(sim41$X7[1:6],sim42$X7[1:6],sim43$X7[1:6]),
                     x = factor(rep(c(150,1200,6000),each=6)))

greeks <- list(bquote(beta[0]), bquote(beta[1]), bquote(beta[2]),bquote(beta[3]),bquote(beta[4]),bquote(phi))

pdf("abserror.pdf")

error1 <- ggplot(data=erdat1, aes(x,y,colour=Coefficient)) + geom_line(aes(group=Coefficient),linetype=3) + geom_point() +
  labs(x="Sample Size",y="Absolute Error") + ggtitle(expression(paste("Asolute Error Plot for ", phi, "=2"))) +
  theme(legend.position="none")
error2 <- ggplot(data=erdat2, aes(x,y,colour=Coefficient)) + geom_line(aes(group=Coefficient),linetype=3) + geom_point() +
  labs(x="Sample Size",y="Absolute Error") + ggtitle(expression(paste("Asolute Error Plot for ", phi, "=5"))) +
  theme(legend.position="none")
error3 <- ggplot(data=erdat3, aes(x,y,colour=Coefficient)) + geom_line(aes(group=Coefficient),linetype=3) + geom_point() +
  labs(x="Sample Size",y="Absolute Error") + ggtitle(expression(paste("Asolute Error Plot for ", phi, "=10"))) +
  theme(legend.position="none")
error4 <- ggplot(data=erdat4, aes(x,y,colour=Coefficient)) + geom_line(aes(group=Coefficient),linetype=3) + geom_point() +
  labs(x="Sample Size",y="Absolute Error") + ggtitle(expression(paste("Asolute Error Plot for ", phi, "=60"))) +
  scale_colour_manual(values=1:6,labels=greeks) + theme(legend.position = c(0.8, 0.6))

grid.arrange(error1,error2,error3,error4,ncol=2)

dev.off()


##########################################################
######### Response Variable Density and QQ Plot ##########
##########################################################

EnvWomen <- read.xls(("environment_data.xlsx"),sheet=1,header=TRUE)

pdf("responsedensity.pdf")

dens <- ggplot(data=EnvWomen,aes(x=EnvDeaths)) + geom_density() + labs(x="Environmental Deaths (%)", main="Density Plot") + xlim(0,0.5)
qq <- ggplot(data=EnvWomen,aes(sample=EnvDeaths)) + geom_qq() + labs(main="QQ Plot")

grid.arrange(dens,qq,ncol=2)
dev.off()

##################################################
######### Comparison of Residuals Plots ##########
##################################################

x= seq(1,length(fit.dat$residuals))
resp1=data_logit$residuals
resp2=weighted.residuals(fit.dat.n)
resp3=weighted.residuals(fit.dat.b)
respbet <- data.frame(x,resp1)
resplin <- data.frame(x,resp2)
respbin <- data.frame(x,resp3)

pdf("residuals.pdf")

resp.beta <- ggplot(respbet) + geom_point(aes(x=x,y=resp1)) +  geom_hline(yintercept = 0) + labs(x="Beta Predictions",y="Residuals") + ylim(-.5,.5) 
resp.lm <- ggplot(resplin) + geom_point(aes(x=x,y=resp2)) +  geom_hline(yintercept = 0) + labs(x="Linear Predictions",y="") + ylim(-.5,.5) 
resp.bin <- ggplot(respbin) + geom_point(aes(x=x,y=resp3)) +  geom_hline(yintercept = 0) + labs(x="Binomial Predictions",y="") + ylim(-.5,.5)

grid.arrange(resp.beta,resp.lm,resp.bin,ncol=3)

dev.off()


