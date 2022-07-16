# Load Packages
library(ri2)
library(ggplot2)
library(xtable)

# Read in Data
ureport <- read.csv(file="ureportmaster.csv",header=T)
View(ureport)

# Treatment=treat_textsearch
# Responses=totalresponses, resp_rate, qual_grade_mean

with(ureport, table(treat_textsearch))
# Control=32533
# HumanFB=345
# ResponseLot=429
# Twitter=380

# Subset Dataset with each Treatment Level with Control
ureportFB <- ureport[which(ureport$treat_textsearch=="Control" | ureport$treat_textsearch=="HumanFB"),]
ureportFB$treat_textsearch <- ifelse(ureportFB$treat_textsearch=="Control",0,1)
ureportFB$treat_textsearch <- as.factor(ureportFB$treat_textsearch)
# N=32878
ureportRL <- ureport[which(ureport$treat_textsearch=="Control" | ureport$treat_textsearch=="ResponseLot"),]
ureportRL$treat_textsearch <- as.factor(ureportRL$treat_textsearch)
ureportRL$treat_textsearch <- ifelse(ureportRL$treat_textsearch=="Control",0,1)
# N=32962
ureportT <- ureport[which(ureport$treat_textsearch=="Control" | ureport$treat_textsearch=="Twitter"),]
ureportT$treat_textsearch <- as.factor(ureportT$treat_textsearch)
ureportT$treat_textsearch <- ifelse(ureportT$treat_textsearch=="Control",0,1)
# N=32913

###################################################
######## Plot the Response with Treatments ########
###################################################

ureport <- ureport[-which(ureport$treat_textsearch==""),]

############################################################
# Response: totalresponses

mean.tr <- aggregate(totalresponses ~ treat_textsearch, data=ureport, mean)
sd.tr <- aggregate(totalresponses ~ treat_textsearch, data=ureport, sd)
n.tr <- aggregate(totalresponses ~ treat_textsearch, data=ureport, length)
dat.tr <- setNames(mean.tr, c("Treatment", "TotalResponses"))

ci.trc <- mean.tr[1,2] + c(-1,1)*qnorm(0.975)*sd.tr[1,2]/sqrt(n.tr[1,2])
ci.trfb <- mean.tr[2,2] + c(-1,1)*qnorm(0.975)*sd.tr[2,2]/sqrt(n.tr[2,2])
ci.trrl <- mean.tr[3,2] + c(-1,1)*qnorm(0.975)*sd.tr[3,2]/sqrt(n.tr[3,2])
ci.trt <- mean.tr[4,2] + c(-1,1)*qnorm(0.975)*sd.tr[4,2]/sqrt(n.tr[4,2])
int.trc <- (ci.trc[2]-ci.trc[1])/2
int.trfb <- (ci.trfb[2]-ci.trfb[1])/2
int.trrl <- (ci.trrl[2]-ci.trrl[1])/2
int.trt <- (ci.trt[2]-ci.trt[1])/2

png("totalresp.png")
ggplot(dat.tr, aes(x=Treatment, y=TotalResponses)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=c(TotalResponses[1]-int.trc,TotalResponses[2]-int.trfb,TotalResponses[3]-int.trrl,TotalResponses[4]-int.trt),
                    ymax=c(TotalResponses[1]+int.trc,TotalResponses[2]+int.trfb,TotalResponses[3]+int.trrl,TotalResponses[4]+int.trt)),
                colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=TotalResponses))
dev.off()

############################################################
# Response: resp_rate

mean.rr <- aggregate(resp_rate ~ treat_textsearch, data=ureport, mean)
sd.rr <- aggregate(resp_rate ~ treat_textsearch, data=ureport, sd)
n.rr <- aggregate(resp_rate ~ treat_textsearch, data=ureport, length)
dat.rr <- setNames(mean.rr, c("Treatment", "ResponseRates"))

ci.rrc <- mean.rr[1,2] + c(-1,1)*qnorm(0.975)*sd.rr[1,2]/sqrt(n.rr[1,2])
ci.rrfb <- mean.rr[2,2] + c(-1,1)*qnorm(0.975)*sd.rr[2,2]/sqrt(n.rr[2,2])
ci.rrrl <- mean.rr[3,2] + c(-1,1)*qnorm(0.975)*sd.rr[3,2]/sqrt(n.rr[3,2])
ci.rrt <- mean.rr[4,2] + c(-1,1)*qnorm(0.975)*sd.rr[4,2]/sqrt(n.rr[4,2])
int.rrc <- (ci.rrc[2]-ci.rrc[1])/2
int.rrfb <- (ci.rrfb[2]-ci.rrfb[1])/2
int.rrrl <- (ci.rrrl[2]-ci.rrrl[1])/2
int.rrt <- (ci.rrt[2]-ci.rrt[1])/2

png("resprate.png")
ggplot(dat.rr, aes(x=Treatment, y=ResponseRates)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=c(ResponseRates[1]-int.rrc,ResponseRates[2]-int.rrfb,ResponseRates[3]-int.rrrl,ResponseRates[4]-int.rrt),
                    ymax=c(ResponseRates[1]+int.rrc,ResponseRates[2]+int.rrfb,ResponseRates[3]+int.rrrl,ResponseRates[4]+int.rrt)),
                colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=ResponseRates))
dev.off()

############################################################
# Response: qual_grade_mean

mean.q <- aggregate(qual_grade_mean ~ treat_textsearch, data=ureport, mean)
sd.q <- aggregate(qual_grade_mean ~ treat_textsearch, data=ureport, sd)
n.q <- aggregate(qual_grade_mean ~ treat_textsearch, data=ureport, length)
dat.q <- setNames(mean.q, c("Treatment", "QualitativeGrade"))

ci.qc <- mean.q[1,2] + c(-1,1)*qnorm(0.975)*sd.q[1,2]/sqrt(n.q[1,2])
ci.qfb <- mean.q[2,2] + c(-1,1)*qnorm(0.975)*sd.q[2,2]/sqrt(n.q[2,2])
ci.qrl <- mean.q[3,2] + c(-1,1)*qnorm(0.975)*sd.q[3,2]/sqrt(n.q[3,2])
ci.qt <- mean.q[4,2] + c(-1,1)*qnorm(0.975)*sd.q[4,2]/sqrt(n.q[4,2])
int.qc <- (ci.qc[2]-ci.qc[1])/2
int.qfb <- (ci.qfb[2]-ci.qfb[1])/2
int.qrl <- (ci.qrl[2]-ci.qrl[1])/2
int.qt <- (ci.qt[2]-ci.qt[1])/2

png("qual.png")
ggplot(dat.q, aes(x=Treatment, y=QualitativeGrade)) +
  geom_bar(position="dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=c(QualitativeGrade[1]-int.qc,QualitativeGrade[2]-int.qfb,QualitativeGrade[3]-int.qrl,QualitativeGrade[4]-int.qt),
                    ymax=c(QualitativeGrade[1]+int.qc,QualitativeGrade[2]+int.qfb,QualitativeGrade[3]+int.qrl,QualitativeGrade[4]+int.qt)),
                colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=QualitativeGrade))
dev.off()


#########################################
######## Randomization Inference ########
#########################################

############################################################
# Response: totalresponses

# Treatment=HumanFB:
with(ureportFB, table(treat_textsearch))
dec <- declare_ra(N = 32878, m = 345)

riFB_totres <- conduct_ri(formula = totalresponses ~ treat_textsearch, declaration = dec, 
                  assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportFB)
summary(riFB_totres) 
png("rifb_totres.png")
plot(riFB_totres)
dev.off()

# Treatment=ResponseLot:
with(ureportRL, table(treat_textsearch))
dec <- declare_ra(N = 32962, m = 429)

riRL_totres <- conduct_ri(formula = totalresponses ~ treat_textsearch, declaration = dec, 
                          assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportRL)
summary(riRL_totres) 
png("rirl_totres.png")
plot(riRL_totres)
dev.off()

# Treatment=Twitter:
with(ureportT, table(treat_textsearch))
dec <- declare_ra(N = 32913, m = 380)

riT_totres <- conduct_ri(formula = totalresponses ~ treat_textsearch, declaration = dec, 
                          assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportT)
summary(riT_totres) 
png("rit_totres.png")
plot(riT_totres)
dev.off()

############################################################
# Response: resp_rate

# Remove observations without reponse rate
ureportr <- ureport[-which(is.na(ureport$resp_rate)),]
# Subset Dataset with each Treatment Level with Control
ureportrFB <- ureportr[which(ureportr$treat_textsearch=="Control" | ureportr$treat_textsearch=="HumanFB"),]
ureportrFB$treat_textsearch <- ifelse(ureportrFB$treat_textsearch=="Control",0,1)
ureportrFB$treat_textsearch <- as.factor(ureportrFB$treat_textsearch)
# N=32877
ureportrRL <- ureportr[which(ureportr$treat_textsearch=="Control" | ureportr$treat_textsearch=="ResponseLot"),]
ureportrRL$treat_textsearch <- as.factor(ureportrRL$treat_textsearch)
ureportrRL$treat_textsearch <- ifelse(ureportrRL$treat_textsearch=="Control",0,1)
# N=32961
ureportrT <- ureportr[which(ureportr$treat_textsearch=="Control" | ureportr$treat_textsearch=="Twitter"),]
ureportrT$treat_textsearch <- as.factor(ureportrT$treat_textsearch)
ureportrT$treat_textsearch <- ifelse(ureportrT$treat_textsearch=="Control",0,1)
# N=32913

# Treatment=HumanFB:
with(ureportFB, table(treat_textsearch))
dec <- declare_ra(N = 32877, m = 345)

riFB_resrate <- conduct_ri(formula = resp_rate ~ treat_textsearch, declaration = dec, 
                          assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportrFB)
summary(riFB_resrate) 
png("rifb_resrate.png")
plot(riFB_resrate)
dev.off()

# Treatment=ResponseLot:
with(ureportRL, table(treat_textsearch))
dec <- declare_ra(N = 32962, m = 429)

riRL_resrate <- conduct_ri(formula = resp_rate ~ treat_textsearch, declaration = dec, 
                          assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportrRL)
summary(riRL_resrate)
png("rirl_resrate.png")
plot(riRL_resrate)
dev.off()

# Treatment=Twitter:
with(ureportT, table(treat_textsearch))
dec <- declare_ra(N = 32913, m = 380)

riT_resrate <- conduct_ri(formula = resp_rate ~ treat_textsearch, declaration = dec, 
                         assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportrT)
summary(riT_resrate) 
png("rit_resrate.png")
plot(riT_resrate)
dev.off()

############################################################
# Response: qual_grade_mean

# Remove observations without qualitative grade
ureportqual <- ureport[-which(is.na(ureport$qual_grade_mean)),]
# Subset Dataset with each Treatment Level with Control
ureportqFB <- ureportqual[which(ureportqual$treat_textsearch=="Control" | ureportqual$treat_textsearch=="HumanFB"),]
ureportqFB$treat_textsearch <- ifelse(ureportqFB$treat_textsearch=="Control",0,1)
ureportqFB$treat_textsearch <- as.factor(ureportqFB$treat_textsearch)
# N=2781
ureportqRL <- ureportqual[which(ureportqual$treat_textsearch=="Control" | ureportqual$treat_textsearch=="ResponseLot"),]
ureportqRL$treat_textsearch <- as.factor(ureportqRL$treat_textsearch)
ureportqRL$treat_textsearch <- ifelse(ureportqRL$treat_textsearch=="Control",0,1)
# N=2858
ureportqT <- ureportqual[which(ureportqual$treat_textsearch=="Control" | ureportqual$treat_textsearch=="Twitter"),]
ureportqT$treat_textsearch <- as.factor(ureportqT$treat_textsearch)
ureportqT$treat_textsearch <- ifelse(ureportqT$treat_textsearch=="Control",0,1)
# N=2819

# Treatment=HumanFB:
with(ureportqFB, table(treat_textsearch))
dec <- declare_ra(N = 2781, m = 76)

riFB_qual <- conduct_ri(formula = qual_grade_mean ~ treat_textsearch, declaration = dec, 
                          assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportqFB)
summary(riFB_qual)
png("rifb_qual.png")
plot(riFB_qual)
dev.off()

# Treatment=ResponseLot:
with(ureportqRL, table(treat_textsearch))
dec <- declare_ra(N = 2858, m = 153)

riRL_qual <- conduct_ri(formula = qual_grade_mean ~ treat_textsearch, declaration = dec, 
                          assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportqRL)
summary(riRL_qual) 
png("rirl_qual.png")
plot(riRL_qual)
dev.off()

# Treatment=Twitter:
with(ureportqT, table(treat_textsearch))
dec <- declare_ra(N = 2819, m = 114)

riT_qual <- conduct_ri(formula = qual_grade_mean ~ treat_textsearch, declaration = dec, 
                         assignment = "treat_textsearch", sharp_hypothesis = 0, sims=10000, data = ureportqT)
summary(riT_qual) 
png("rit_qual.png")
plot(riT_qual)
dev.off()

#####################################
######## Regression Analysis ########
#####################################

ureport$orangetower <- as.factor(ureport$orangetower)
ureport$treat_textsearch <- as.factor(ureport$treat_textsearch)
ureport$ureporter100 <- as.factor(ifelse(ureport$ureporters>=100,1,0))

# sdindex, ureporters, orangetower, 

modelTR <- lm(totalresponses ~ -1 + treat_textsearch + orangetower + ureporters + ureporter100, data=ureport)
summary(modelTR)

xtable(summary(modelTR)$coef)

modelRR <- lm(resp_rate ~ -1 + treat_textsearch + orangetower + ureporters + ureporter100, data=ureport)
summary(modelRR)

xtable(summary(modelRR)$coef)

modelQ <- lm(qual_grade_mean ~ -1 + treat_textsearch + orangetower + ureporters + ureporter100, data=ureport)
summary(modelQ)

xtable(summary(modelQ)$coef)




