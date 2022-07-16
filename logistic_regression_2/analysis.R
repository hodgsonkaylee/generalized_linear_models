####################################################
######### STAT 536  - Tulips in Netherlands ########
####################################################

library(ggplot2)
library(dplyr)
library(gridExtra)
library(nlme)
library(splines)
library(xtable)

#########################################################
#########################################################

# Import Data
tulips <- read.csv(file="tulips_data.csv",header=T)
tulips$Germinated <- as.factor(tulips$Germinated)
tulips$germinated <- ifelse(tulips$Germinated=="N",0,1)
germinated <- ifelse(tulips$Germinated=="N",0,1)
tulips$Population <- as.factor(tulips$Population)
tulips$chillingtime <- tulips$ChillingTime
tulips$ChillingTime <- as.factor(tulips$ChillingTime)
tulipsfin <- tulips[-which(tulips$Population==12),]
View(tulipsfin)

# Subset the data for each population
tulips1 <- tulips[which(tulips$Population==1),]
tulips2 <- tulips[which(tulips$Population==2),]
tulips3 <- tulips[which(tulips$Population==3),]
tulips4 <- tulips[which(tulips$Population==4),]
tulips5 <- tulips[which(tulips$Population==5),]
tulips6 <- tulips[which(tulips$Population==6),]
tulips7 <- tulips[which(tulips$Population==7),]
tulips8 <- tulips[which(tulips$Population==8),]
tulips9 <- tulips[which(tulips$Population==9),]
tulips10 <- tulips[which(tulips$Population==10),]
tulips11 <- tulips[which(tulips$Population==11),]
tulips12 <- tulips[which(tulips$Population==12),]

#########################################################
#########################################################

# Explore Data

# Proportion that germinated for each population
ggplot(tulipsfin) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
plot(tulips$chillingtime,tulips$Germinated)
ggplot(tulips) + geom_point(aes(x = ChillingTime, y= Germinated)) + geom_smooth(aes(x = as.numeric(ChillingTime), y= as.numeric(Germinated)), method = 'loess', se = FALSE)
ggplot(tulips) + geom_point(aes(x = Population, y= Germinated)) + geom_smooth(aes(x = as.numeric(Population), y= as.numeric(Germinated)), se = FALSE)
scatter.smooth(as.numeric(tulips$ChillingTime),tulips$Germinated,degree=2)

# Proportion of germinated for different chilling times for individual populations
ggplot(tulips1) + geom_point(aes(x = ChillingTime, y= Germinated)) + geom_smooth(aes(x = as.numeric(ChillingTime), y= as.numeric(Germinated)), method = 'loess', se = FALSE)
ggplot(tulips2) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips3) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips4) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips5) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips6) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips7) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips8) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips9) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips10) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips11) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")
ggplot(tulips12) + geom_bar(aes(x = ChillingTime, fill = Germinated), width=.5, position = "dodge")

# General Overall Trends 
pdf("GermChill.pdf")
par(mfrow=c(1,2))
plot(tulips$ChillingTime,tulips$Germinated,
     ylab="Germinated",xlab="Chilling Time",
     col=c("darkmagenta","cyan3"))
scatter.smooth(as.numeric(tulips$ChillingTime),tulips$Germinated,degree=2,
                ylab="Germinated",xlab="Chilling Time",yaxt = 'n',
               col=ifelse(tulips$Germinated=="N","darkmagenta","cyan3"))
axis(side=2,at=c(1,2),labels=c('No','Yes'),las=2)
dev.off()

xtabs(formula=~ChillingTime+Germinated, data=tulips)

tab1<- table(tulips$ChillingTime,tulips$Germinated)
ftable(prop.table(tab1,1))

prop.table(tulips$ChillingTime+tulips$Germinated)

# Chilling time for each population at each different chilling time
descript <- matrix(NA,nrow=12,ncol=7)
for(i in 1:7){
  descript[1,i] <- mean(tulips1$germinated[which(tulips1$ChillingTime==(i-1)*2)])
  descript[2,i] <- mean(tulips2$germinated[which(tulips2$ChillingTime==(i-1)*2)])
  descript[3,i] <- mean(tulips3$germinated[which(tulips3$ChillingTime==(i-1)*2)])
  descript[4,i] <- mean(tulips4$germinated[which(tulips4$ChillingTime==(i-1)*2)])
  descript[5,i] <- mean(tulips5$germinated[which(tulips5$ChillingTime==(i-1)*2)])
  descript[6,i] <- mean(tulips6$germinated[which(tulips6$ChillingTime==(i-1)*2)])
  descript[7,i] <- mean(tulips7$germinated[which(tulips7$ChillingTime==(i-1)*2)])
  descript[8,i] <- mean(tulips8$germinated[which(tulips8$ChillingTime==(i-1)*2)])
  descript[9,i] <- mean(tulips9$germinated[which(tulips9$ChillingTime==(i-1)*2)])
  descript[10,i] <- mean(tulips10$germinated[which(tulips10$ChillingTime==(i-1)*2)])
  descript[11,i] <- mean(tulips11$germinated[which(tulips11$ChillingTime==(i-1)*2)])
  descript[12,i] <- mean(tulips12$germinated[which(tulips12$ChillingTime==(i-1)*2)])
}

xtable(descript)

#########################################################
#########################################################
# Fit with ChillingTime as continuous, with splines, vs. Fit with ChillingTime as Categorical
# Cross-validation to find best number of knots for natural splines (0:5) or factor
accuracy <- matrix(NA,nrow=1000,ncol=7)
accuracyos <- matrix(NA,nrow=1000,ncol=7)
bestacc <- rep(NA,1000)
bestaccos <- rep(NA,1000)
for(i in 1:1000){
  x <- sample(dim(tulipsfin)[1],252)
  test <- tulipsfin[x,]
  train <- tulipsfin[-x,]
  for(j in 1:6){
    testgls <- glm(germinated ~ -1 + ns(chillingtime,df=j):Population, family=binomial(link="logit"), data=train)
    # Overall Accuracy In Sample
    pred.prob <- predict(testgls,type="response")
    pred.class <- rep(0,length(pred.prob))
    pred.class[pred.prob>0.5] <- 1
    accuracy[i,j] <- mean(pred.class==train$germinated)
    # Overall Accuracy Out of Sample
    pred.probos <- predict(testgls,test,type="response")
    pred.classos <- rep(0,length(pred.probos))
    pred.classos[pred.probos>0.5] <- 1
    accuracyos[i,j] <- mean(pred.classos==test$germinated)
  }
  testglsf <- glm(Germinated ~ -1 + ChillingTime:Population, family=binomial(link="logit"), data=train)
  # Overall Accuracy In Sample
  pred.prob <- predict(testglsf,type="response")
  pred.class <- rep(0,length(pred.prob))
  pred.class[pred.prob>0.5] <- 1
  accuracy[i,7] <- mean(pred.class==train$germinated)
  # Best Accuracy In Sample
  bestacc[i] <- which.max(accuracy[,j])
  # Overall Accuracy Out of Sample
  pred.probos <- predict(testglsf,test,type="response")
  pred.classos <- rep(0,length(pred.probos))
  pred.classos[pred.probos>0.5] <- 1
  accuracyos[i,7] <- mean(pred.classos==test$germinated)
}

# Estimate the number of Times each model had the best/tied best overall accuracy
bestacc <- matrix(NA,nrow=1000,ncol=7)
for(i in 1:1000){
  for(j in 1:7){
  bestacc[i,j] <- ifelse(max(accuracy[i,])==accuracy[i,j],TRUE,FALSE)
  }
}
bestaccos <- matrix(NA,nrow=1000,ncol=7)
for(i in 1:1000){
  for(j in 1:7){
    bestaccos[i,j] <- ifelse(max(accuracyos[i,])==accuracyos[i,j],TRUE,FALSE)
  }
}
bestaccuracy <- apply(bestacc,2,sum)
bestaccuracyos <- apply(bestaccos,2,sum)
bestac <- c(rep(2,0),rep(3,1),rep(4,3),rep(5,556),rep(6,990),rep(7,1000))
apply(accuracy,2,mean)
apply(accuracyos,2,mean)

pdf("tuning.pdf")
par(mfrow=c(1,2))
plot(x=1:7,y=apply(accuracy,2,mean),main="Average Accuracy for Knots=1:5",xlab="Knots",ylab="Average Overall Accuracy",xaxt='n')
abline(v=7,col='blue')
axis(1, at = 1:7, label = c(1:6,"Factor"), las=2, tick=F, cex.axis=1.25)
hist(bestac,breaks=1:7,main="Frequency of Optimal Values",xlab="Knots",xaxt='n')
axis(1, at = 1:7, label = c(1:6,"Factor"), las=2, tick=F, cex.axis=1.25)
dev.off()

##############################################################
# Fit with ChillingTime as a factor
fit1 <- glm(Germinated ~ -1 + ChillingTime:Population, family=binomial(link="logit"), data=tulipsfin)
summary(fit1)

pred.prob <- predict(fit1,type="response")
pred.class <- rep(0,length(pred.prob))
pred.class[pred.prob>0.5] <- 1
mean(pred.class==tulipsfin$germinated)

library(ROCR)
pred <- prediction(predictions=pred.prob,labels=tulipsfin$germinated)
rocr.perf <-performance(pred,'tpr','fpr')
pdf("ROCcurve.pdf")
plot(rocr.perf)
abline(a=0,b=1,ylim=c(0,1))
dev.off()

rocr.acur <-performance(pred,'acc')
max.acur <- unlist(rocr.acur@x.values)[which.max(unlist(rocr.acur@y.values))]
pdf("accuracyplot.pdf")
plot(rocr.acur)
abline(v=max.acur)
dev.off()
# Cut-off value = 0.5177266

# Model Fit and Prediction Power
pred <- prediction(predictions=pred.prob,labels=tulipsfin$germinated)
performance(pred,measure='auc') # AUC=0.8829519
unlist(performance(pred,measure='spec')@y.values)[which(round(unlist(performance(pred,measure='spec')@x.values),2)==0.52)[2]] # Specificity=0.77
unlist(performance(pred,measure='sens')@y.values)[which(round(unlist(performance(pred,measure='sens')@x.values),2)==0.52)[2]] # Sensitivity=0.83
unlist(performance(pred,measure='ppv')@y.values)[which(round(unlist(performance(pred,measure='ppv')@x.values),2)==0.52)[2]]
unlist(performance(pred,measure='npv')@y.values)[which(round(unlist(performance(pred,measure='npv')@x.values),2)==0.52)[2]]
unlist(performance(pred,measure='acc')@y.values)[which(round(unlist(performance(pred,measure='acc')@x.values),2)==0.52)[2]]

modfitcv <- matrix(NA,nrow=1000,ncol=6)
for(i in 1:1000){
  x <- sample(dim(tulipsfin)[1],252)
  test <- tulipsfin[x,]
  train <- tulipsfin[-x,]
  testgls <- glm(germinated ~ -1 + ChillingTime:Population, family=binomial(link="logit"), data=train)
  pred.prob <- predict(testgls,test,type="response")
  pred <- prediction(predictions=pred.prob,labels=test$germinated)
  rocr.acur <-performance(pred,'acc')
  max.acur <- unlist(rocr.acur@x.values)[which.max(unlist(rocr.acur@y.values))]
  modfitcv[i,1] <- unlist(performance(pred,measure='auc')@y.values)
  modfitcv[i,2] <- unlist(performance(pred,measure='spec')@y.values)[which(round(unlist(performance(pred,measure='spec')@x.values),2)==round(max.acur,2))[2]]
  modfitcv[i,3] <- unlist(performance(pred,measure='sens')@y.values)[which(round(unlist(performance(pred,measure='sens')@x.values),2)==round(max.acur,2))[2]]
  modfitcv[i,4] <- unlist(performance(pred,measure='ppv')@y.values)[which(round(unlist(performance(pred,measure='ppv')@x.values),2)==round(max.acur,2))[2]]
  modfitcv[i,5] <- unlist(performance(pred,measure='npv')@y.values)[which(round(unlist(performance(pred,measure='npv')@x.values),2)==round(max.acur,2))[2]]
  modfitcv[i,6] <- unlist(performance(pred,measure='acc')@y.values)[which(round(unlist(performance(pred,measure='acc')@x.values),2)==round(max.acur,2))[2]]
}
apply(modfitcv,2,mean,na.rm=TRUE) # CV Values: 0.8644369 0.8247139 0.7678527 0.7667072 0.8331383 0.7958162

#################################################################
# Model Results
fit1 <- glm(Germinated ~ -1 + ChillingTime:Population, family=binomial(link="logit"), data=tulipsfin)
sum.glm <- summary(fit1)
upperci <- sum.glm$coef[,1] + qnorm(0.997)*sum.glm$coef[,2] 
lowerci <- sum.glm$coef[,1] - qnorm(0.997)*sum.glm$coef[,2]
xtable(cbind(sum.glm$coef[1:14,1],lowerci[1:14],upperci[1:14],sum.glm$coef[1:14,4]))

#################################################################
# Question 1:Is the probability of germination for each chilling time the same across all populations? Which populations are same/different?

fit1r <- glm(Germinated ~ -1 + ChillingTime, family=binomial(link="logit"), data=tulipsfin)
anova(fit1r,fit1,test="LRT")
pred.prob <- predict(fit1r,type="response")
pred.class <- rep(0,length(pred.prob))
pred.class[pred.prob>0.5] <- 1

# Compare Probabilities
# Predicted Probabilities:
summary(fit1)
pred.prob <- predict(fit1,type="response")
pred.class <- rep(0,length(pred.prob))
pred.class[pred.prob>0.5] <- 1

options(scipen=999)
probs <- matrix(NA,nrow=12,ncol=7)
for(i in 1:11){
  probs[i,] <- c(mean(pred.prob[which(tulipsfin$Population==i & tulipsfin$ChillingTime==0)]),
                 mean(pred.prob[which(tulipsfin$Population==i & tulipsfin$ChillingTime==2)]),
                 mean(pred.prob[which(tulipsfin$Population==i & tulipsfin$ChillingTime==4)]),
                 mean(pred.prob[which(tulipsfin$Population==i & tulipsfin$ChillingTime==6)]),
                 mean(pred.prob[which(tulipsfin$Population==i & tulipsfin$ChillingTime==8)]),
                 mean(pred.prob[which(tulipsfin$Population==i & tulipsfin$ChillingTime==10)]),
                 mean(pred.prob[which(tulipsfin$Population==i & tulipsfin$ChillingTime==12)]))
}
probs[12,] <- cbind(mean(pred.prob[tulipsfin$ChillingTime==0]),
      mean(pred.prob[tulipsfin$ChillingTime==2]),
      mean(pred.prob[tulipsfin$ChillingTime==4]),
      mean(pred.prob[tulipsfin$ChillingTime==6]),
      mean(pred.prob[tulipsfin$ChillingTime==8]),
      mean(pred.prob[tulipsfin$ChillingTime==10]),
      mean(pred.prob[tulipsfin$ChillingTime==12]))

xtable(probs)

pdf("q1plots.pdf")
par(mfrow=c(3,4))
plot(1:7,probs[1,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 1",xaxt='n',col='blue')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[2,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 2",xaxt='n',col='blue')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[3,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 3",xaxt='n',col='blue')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[4,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 4",xaxt='n',col='purple')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[5,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 5",xaxt='n',col='red')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[6,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 6",xaxt='n',col='purple')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[7,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 7",xaxt='n',col='purple')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[8,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 8",xaxt='n',col='blue')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[9,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 9",xaxt='n',col='darkgreen')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[10,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 10",xaxt='n',col='darkgreen')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[11,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 11",xaxt='n',col='purple')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
dev.off()

pdf("q2plot.pdf")
plot(1:7,probs[12,],ylab="Probability of Germination",xlab="Chilling Time",xaxt='n')
abline(v=6,col='blue')
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=1.5)
dev.off()

pdf("q2aplots.pdf")
par(mfrow=c(3,4))
plot(1:7,probs[1,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 1",xaxt='n',col='blue')
abline(v=6)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[2,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 2",xaxt='n',col='blue')
abline(v=6)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[3,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 3",xaxt='n',col='purple')
abline(v=5)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[4,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 4",xaxt='n',col='blue')
abline(v=6)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[5,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 5",xaxt='n',col='red')
abline(v=2)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[6,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 6",xaxt='n',col='blue')
abline(v=6)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[7,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 7",xaxt='n',col='blue')
abline(v=6)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[8,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 8",xaxt='n',col='purple')
abline(v=5)
abline(v=4)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[9,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 9",xaxt='n',col='blue')
abline(v=6)
abline(v=7)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[10,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 10",xaxt='n',col='blue')
abline(v=6)
abline(v=5)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
plot(1:7,probs[11,],ylab="Probability of Germination",xlab="Chilling Time",main="Population 11",xaxt='n',col='blue')
abline(v=6)
axis(1, at = 1:7, label = (0:6)*2, las=1, cex.axis=.5)
dev.off()

# F test for 3rd question:
C <- rbind(rep(c(0,0,0,0,-1/11,1/11,0),11))
bhat <- coef(fit1)
df1 <- 2310-78
df2 <- 1
vhatbhat <- vcov(fit1)
F.val <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val
p.val <- 1-pf(F.val,df2,df1)
p.val

pdf("q3plot.pdf")
plot(x=1:2,y=c(2,1),type="b",xlab="Chilling Time (weeks)",ylab="Probability of Germination",xaxt='n',yaxt='n',col="purple")
axis(1, at = 1:2, label = c(10,8), las=1, cex.axis=1.5)
axis(2, at=1:2, label=c("63%","77%"), las=1, cex.axis=1.5)
dev.off()

# F test for 3rd question - individual populations:
C <- rbind(c(c(0,0,0,0,-1,1,0),rep(0,70)))
F.val1 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val1
p.val1 <- 1-pf(F.val1,df2,df1)
p.val1

C <- rbind(c(rep(0,7),c(0,0,0,0,-1,1,0),rep(0,63)))
F.val2 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val2
p.val2 <- 1-pf(F.val2,df2,df1)
p.val2

C <- rbind(c(rep(0,14),c(0,0,0,0,-1,1,0),rep(0,56)))
F.val3 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val3
p.val3 <- 1-pf(F.val3,df2,df1)
p.val3

C <- rbind(c(rep(0,21),c(0,0,0,0,-1,1,0),rep(0,49)))
F.val4 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val4
p.val4 <- 1-pf(F.val4,df2,df1)
p.val4

C <- rbind(c(rep(0,28),c(0,0,0,0,-1,1,0),rep(0,42)))
F.val5 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val5
p.val5 <- 1-pf(F.val5,df2,df1)
p.val5

C <- rbind(c(rep(0,35),c(0,0,0,0,-1,1,0),rep(0,35)))
F.val6 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val6
p.val6 <- 1-pf(F.val6,df2,df1)
p.val6

C <- rbind(c(rep(0,42),c(0,0,0,0,-1,1,0),rep(0,28)))
F.val7 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val7
p.val7 <- 1-pf(F.val7,df2,df1)
p.val7

C <- rbind(c(rep(0,49),c(0,0,0,0,-1,1,0),rep(0,21)))
F.val8 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val8
p.val8 <- 1-pf(F.val8,df2,df1)
p.val8

C <- rbind(c(rep(0,56),c(0,0,0,0,-1,1,0),rep(0,14)))
F.val9 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val9
p.val9 <- 1-pf(F.val9,df2,df1)
p.val9

C <- rbind(c(rep(0,63),c(0,0,0,0,-1,1,0),rep(0,7)))
F.val10 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val10
p.val10 <- 1-pf(F.val10,df2,df1)
p.val10

C <- rbind(c(rep(0,70),c(0,0,0,0,-1,1,0)))
F.val11 <- t(C %*% bhat) %*% solve(C %*% vhatbhat %*% t(C)) %*% C %*% bhat / df2
F.val11
p.val11 <- 1-pf(F.val11,df2,df1)
p.val11

xtable(matrix(c(F.val1,p.val1,
                F.val2,p.val2,
                F.val3,p.val3,
                F.val4,p.val4,
                F.val5,p.val5,
                F.val6,p.val6,
                F.val7,p.val7,
                F.val8,p.val8,
                F.val9,p.val9,
                F.val10,p.val10,
                F.val11,p.val11),byrow=T,nrow=11,ncol=2))

pdf("q3aplot.pdf")
par(mfrow=c(2,2))
plot(x=1:2,y=c(90,97),type="l",xlab="Chilling Time (weeks)",ylab="Probability of Germination",main="Increased",ylim=c(0,100),lty=2,lwd=2,xaxt='n',yaxt='n',col="cadetblue") #3
lines(c(57,70),col='cadetblue1',lty=2,lwd=2) #5
lines(c(30,33),col='cadetblue3',lty=2,lwd=2) #8
axis(1, at = 1:2, label = c(10,8), las=1, cex.axis=1)
axis(2, at=seq(from=0,to=100,by=20), label=c("0%","20%","40%","60%","80%","100%"), las=1, cex.axis=1)
legend('bottomleft',c("3","5","8"),title="Population",col=c("cadetblue","cadetblue1","cadetblue3"),lty=2,lwd=2)
plot(x=1:2,y=c(87,87),type="l",xlab="Chilling Time (weeks)",ylab="Probability of Germination",main="Same",ylim=c(0,100),lty=2,lwd=2,xaxt='n',yaxt='n',col="firebrick") #6
axis(1, at = 1:2, label = c(10,8), las=1, cex.axis=1)
axis(2, at=seq(from=0,to=100,by=20), label=c("0%","20%","40%","60%","80%","100%"), las=1, cex.axis=1)
legend('bottomleft',c("10"),title="Population",col="firebrick",lty=2,lwd=2)
plot(x=1:2,y=c(97,87),type="l",xlab="Chilling Time (weeks)",ylab="Probability of Germination",main="Unsignificantly Decreased",ylim=c(0,100),lty=2,lwd=2,xaxt='n',yaxt='n',col="darkorchid") #1
lines(c(90,83),col='darkorchid1',lty=2,lwd=2) #2
lines(c(90,73),col='darkorchid4',lty=2,lwd=2) #4
lines(c(83,67),col='magenta',lty=2,lwd=2) #11
axis(1, at = 1:2, label = c(10,8), las=1, cex.axis=1)
axis(2, at=seq(from=0,to=100,by=20), label=c("0%","20%","40%","60%","80%","100%"), las=1, cex.axis=1)
legend('bottomleft',c("1","2","4","11"),title="Population",col=c("darkorchid","darkorchid1","darkorchid4","magenta"),lty=2,lwd=2)
plot(x=1:2,y=c(80,43),type="l",xlab="Chilling Time (weeks)",ylab="Probability of Germination",main="Significantly Decreased",ylim=c(0,100),lty=2,lwd=2,xaxt='n',yaxt='n',col="darkolivegreen") #6
lines(c(83,47),col='darkolivegreen3',lty=2,lwd=2) #7
lines(c(60,7),col='darkolivegreen4',lty=2,lwd=2) #9
axis(1, at = 1:2, label = c(10,8), las=1, cex.axis=1)
axis(2, at=seq(from=0,to=100,by=20), label=c("0%","20%","40%","60%","80%","100%"), las=1, cex.axis=1)
legend('bottomleft',c("6","7","9"),title="Population",col=c("darkolivegreen","darkolivegreen3","darkolivegreen4"),lty=2,lwd=2)
dev.off()


