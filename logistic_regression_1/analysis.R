#####################################
####### Import and Clean Data #######
#####################################

iegender <- read.csv('iegenderdat.csv',header=TRUE)
View(iegender)

# Grouped gender of enumerator and gender of official
# 00 = Male enumerator, Male official
# 01 = Male enumerator, Female official
# 10 = Female enumerator, Male official
# 11 = Female enumerator, Female official
genderint <- paste(ifelse(is.na(iegender$enumfemale),"",iegender$enumfemale),ifelse(is.na(iegender$officialfemale), "", iegender$officialfemale),sep='') 
genderint <- ifelse(genderint=="0",NA,genderint)
genderint <- ifelse(genderint=="1",NA,genderint)
iegender$genderint <- as.factor(genderint)

# iegender$socialproof <- ifelse(genderint=="-",NA,iegender$socialproof)
iegender$individualpriority <- ifelse(iegender$individualpriority=="",NA,iegender$individualpriority)
iegender$individualpriority <- ifelse(iegender$individualpriority=="\\\\",NA,iegender$individualpriority)


# Dataframes for each country #
iegenderIndia <- iegender[1:9239,]
iegenderTanz <- iegender[9240:10549,]
iegenderPeru <- iegender[10550:14086,]

# Indicators for Country in overall dataset #
iegender$Country <- as.factor(iegender$Country)
iegender$Indiai <- model.matrix(~-1 + Country, data=iegender)[,1]
iegender$Tanzi <- model.matrix(~-1 + Country, data=iegender)[,2]
iegender$Perui <- model.matrix(~-1 + Country, data=iegender)[,3]

# Group enumerator gender by country
with(iegender, table(Country, enumfemale))
with(iegender, table(Country, socialproof))
# Female Enumerators: 
# India - 5769/9239 (62.44%) 
# Peru - 1496/3537 (42.24%) 
# Tanzania - 658/1310 (50.23%)
# Total - 7923/14086 (56.25%)

# Load Packages # 
library(nlme)
library(lme4)
library(ri2)
library(ggplot2)
library(Rmisc)
library(coefplot)

###################################
####### Difference of Means #######
###################################

# Cumulative difference of means - enumerator gender
t.test(posresponse ~ socialproof, data=iegender) # women performed better, but not significantly better
t.test(held ~ enumfemale, data=iegender) # men performed better, but not significantly better

# Cumulative difference of means - enumerator and official gender
summary(aov(posresponse ~ genderint, data=iegender)) # no significant difference, but interesting gender trend

meanci.ie <- group.CI(posresponse ~ genderint, data=iegender, ci = 0.95)
for(i in 1:4) meanci.ie[i,2] <- meanci.ie[i,2]-meanci.ie[i,3]
dat.ie <- setNames(meanci.ie, c("gender","u", "pos","l"))
pdf("ovrlgendint.pdf")
ggplot(dat.ie, aes(x=gender, y=pos)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin=pos+u,ymax=pos-u),colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=pos)) + theme_light()
dev.off()

summary(aov(held ~ genderint, data=iegender)) # no significant difference - different trends

meanci.ieh <- group.CI(held ~ genderint, data=iegender, ci = 0.95)
for(i in 1:4) meanci.ieh[i,2] <- meanci.ieh[i,2]-meanci.ieh[i,3]
dat.ieh <- setNames(meanci.ieh, c("gender","u", "held","l"))
ggplot(dat.ieh, aes(x=gender, y=held)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin=held+u,ymax=held-u),colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=held)) + theme_light()


# By country difference of means - enumerator gender
## India ##
t.test(posresponse ~ enumfemale, data=iegenderIndia) # not significance - women performed slightly better
t.test(held ~ enumfemale, data=iegenderIndia) # not significant - men performed slightly better
## Tanzania ##
t.test(posresponse ~ enumfemale, data=iegenderTanz) # not significant - women performed better
t.test(held ~ enumfemale, data=iegenderTanz) # not significant - women performed better
## Peru ##
t.test(posresponse ~ enumfemale, data=iegenderPeru) # significant - women performed better
t.test(held ~ enumfemale, data=iegenderPeru) # significant - women performed better

# By country difference of means - enumerator and official gender
## India ##
summary(aov(posresponse ~ genderint, data=iegenderIndia)) # not significant
summary(aov(held ~ genderint, data=iegenderIndia)) # not significant
## Tanzania ##
summary(aov(posresponse ~ genderint, data=iegenderTanz)) # not significant
summary(aov(held ~ genderint, data=iegenderTanz)) # not significant
## Peru ##
summary(aov(posresponse ~ genderint, data=iegenderPeru)) # significant
TukeyHSD(aov(posresponse ~ genderint, data=iegenderPeru)) # women messengers got more pos responses with women than men got with women
summary(aov(held ~ genderint, data=iegenderPeru)) # not significant

meanci.iePp <- group.CI(posresponse ~ enumfemale, data=iegenderPeru, ci = 0.95)
for(i in 1:2) meanci.iePp[i,2] <- meanci.iePp[i,2]-meanci.iePp[i,3]
dat.iePp <- setNames(meanci.iePp, c("gender","u", "pos","l"))

meanci.iePh <- group.CI(held ~ enumfemale, data=iegenderPeru, ci = 0.95)
for(i in 1:2) meanci.iePh[i,2] <- meanci.iePh[i,2]-meanci.iePh[i,3]
dat.iePh <- setNames(meanci.iePh, c("gender","u", "pos","l"))

meanci.ieP <- group.CI(posresponse ~ genderint, data=iegenderPeru, ci = 0.95)
for(i in 1:4) meanci.ieP[i,2] <- meanci.ieP[i,2]-meanci.ieP[i,3]
dat.ieP <- setNames(meanci.ieP, c("gender","u", "pos","l"))

pdf("Perudm.pdf")
par(mfrow=c(1,3))
ggplot(dat.iePp, aes(x=gender, y=pos)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin=pos+u,ymax=pos-u),colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=pos)) + theme_light()
ggplot(dat.iePh, aes(x=gender, y=pos)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin=pos+u,ymax=pos-u),colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=pos)) + theme_light()
ggplot(dat.ieP, aes(x=gender, y=pos)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin=pos+u,ymax=pos-u),colour="black",width = .5) +
  geom_point(position=position_dodge(.9), aes(y=pos)) + theme_light()
dev.off()

#############################################################
####### Difference of Means - Randomization Inference #######
#############################################################

# Cumulative Difference of Means with Randomization Inference
dec <- declare_ra(N = 14086, m = 7923)

ri_posresponse <- conduct_ri(formula = posresponse ~ enumfemale, declaration = dec, 
                      assignment = "enumfemale", sharp_hypothesis = 0, data = iegender)
summary(ri_posresponse) # p-val=.139 - not significant
plot(ri_posresponse)

ri_held <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                      assignment = "enumfemale", sharp_hypothesis = 0, data = iegender)
summary(ri_held) #p-val=.713 - not significant
plot(ri_held)

# By country Difference of Means with Randomization Inference
## India ##
dec <- declare_ra(N = 9239, m = 5769)

ri_posresponseIndia <- conduct_ri(formula = posresponse ~ enumfemale, declaration = dec, 
                             assignment = "enumfemale", sharp_hypothesis = 0, data = iegenderIndia)
summary(ri_posresponseIndia) # p-val=.529 - not significant
plot(ri_posresponseIndia)

ri_heldIndia <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                      assignment = "enumfemale", sharp_hypothesis = 0, data = iegenderIndia)
summary(ri_heldIndia) # p-val=.12 - not significant
plot(ri_heldIndia)

## Tanzania ##
dec <- declare_ra(N = 1310, m = 658)

ri_posresponseTanz <- conduct_ri(formula = posresponse ~ enumfemale, declaration = dec, 
                                  assignment = "enumfemale", sharp_hypothesis = 0, data = iegenderTanz)
summary(ri_posresponseTanz) # p-val=.199 - not significant
plot(ri_posresponseTanz)

ri_heldTanz <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                           assignment = "enumfemale", sharp_hypothesis = 0, data = iegenderTanz)
summary(ri_heldTanz) # p-val=.09 - significant at .1 level
plot(ri_heldTanz)

## Peru ##
dec <- declare_ra(N = 3537, m = 1496)

ri_posresponsePeru <- conduct_ri(formula = posresponse ~ enumfemale, declaration = dec, 
                                 assignment = "enumfemale", sharp_hypothesis = 0, data = iegenderPeru)
summary(ri_posresponsePeru) # p-val=.003
plot(ri_posresponsePeru)

ri_heldPeru <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                          assignment = "enumfemale", sharp_hypothesis = 0, data = iegenderPeru)
summary(ri_heldPeru) # p-val=.027
plot(ri_heldPeru)

### Plot all of the RI distributions ###

pdf("ri_ateov.pdf")
plot(ri_posresponse,main="Overall")
dev.off()
pdf("ri_ateov1.pdf")
plot(ri_held,main="Overall")
dev.off()

pdf("ri_ateP.pdf")
plot(ri_posresponsePeru,main="Peru")
dev.off()
pdf("ri_ateP1.pdf")
plot(ri_heldPeru,main="Peru")
dev.off()

pdf("ri_ate.pdf")
par(mfrow=c(4,2))
plot(ri_posresponse,main="Overall")
plot(ri_held,main="Overall")
plot(ri_posresponseIndia,main="India")
plot(ri_heldIndia,main="India")
plot(ri_posresponseTanz,main="Tanzania")
plot(ri_heldTanz,main="Tanzania")
plot(ri_posresponsePeru,main="Peru")
plot(ri_heldPeru,main="Peru")
dev.off()

###############################################
####### ANOVA - Randomization Inference #######
###############################################

with(iegender, table(Country, genderint))
# Female Enumerators: 
# India - 00:2463, 01:823, 10:4104, 11:1319 
# Peru - 00:1349, 01:673, 10:1003, 11:484
# Tanzania - 00:345, 01:284, 10:341, 11:291
# Total - 00:4157, 01:1780, 10:5448, 11:2094

# Cumulative ANOVA with Randomization Inference
row.has.na <- is.na(iegender$genderint)
iegendred <- iegender[!row.has.na,]

dec <- declare_ra(N = 13479, num_arms = 4, m_each=c(4157,1780,5448,2094),
                  condition_names=c('00','01','10','11'))
anova_pos <-conduct_ri(model_1 = posresponse ~ 1, # restricted model
                      model_2 = posresponse ~ genderint, # unrestricted model
                      assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendred)
summary(anova_pos) # p-value=.293
plot(anova_pos)

anova_held <-conduct_ri(model_1 = held ~ 1, # restricted model
                       model_2 = held ~ genderint, # unrestricted model
                       assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendred)
summary(anova_held) # p-value=.355
plot(anova_held)

# By country ANOVA with Randomization Inference
## India ##
row.has.na <- is.na(iegenderIndia$genderint)
iegendredIndia <- iegenderIndia[!row.has.na,]

dec <- declare_ra(N = 8709, num_arms = 4, m_each=c(2463,823,4104,1319),
                  condition_names=c('00','01','10','11'))
anova_posIndia <-conduct_ri(model_1 = posresponse ~ 1, # restricted model
                       model_2 = posresponse ~ genderint, # unrestricted model
                       assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendredIndia)
summary(anova_posIndia) # p-value=.816
plot(anova_posIndia)

anova_heldIndia <-conduct_ri(model_1 = held ~ 1, # restricted model
                            model_2 = held ~ genderint, # unrestricted model
                            assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendredIndia)
summary(anova_heldIndia) # p-value=.123
plot(anova_heldIndia)

## Peru ##
row.has.na <- is.na(iegenderPeru$genderint)
iegendredPeru <- iegenderPeru[!row.has.na,]

dec <- declare_ra(N = 3509, num_arms = 4, m_each=c(1349,673,1003,484),
                  condition_names=c('00','01','10','11'))
anova_posPeru <-conduct_ri(model_1 = posresponse ~ 1, # restricted model
                            model_2 = posresponse ~ genderint, # unrestricted model
                            assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendredPeru)
summary(anova_posPeru) # p-value=.017
pdf("anovaposperu.pdf")
plot(anova_posPeru)
dev.off()

anova_heldPeru <-conduct_ri(model_1 = held ~ 1, # restricted model
                             model_2 = held ~ genderint, # unrestricted model
                             assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendredPeru)
summary(anova_heldPeru) # p-value=.128
plot(anova_heldPeru)

## Tanzania ##
row.has.na <- is.na(iegenderTanz$genderint)
iegendredTanz <- iegenderTanz[!row.has.na,]

dec <- declare_ra(N = 1261, num_arms = 4, m_each=c(345,284,341,291),
                  condition_names=c('00','01','10','11'))
anova_posTanz <-conduct_ri(model_1 = posresponse ~ 1, # restricted model
                           model_2 = posresponse ~ genderint, # unrestricted model
                           assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendredTanz)
summary(anova_posTanz) # p-value=.338
plot(anova_posTanz)

anova_heldTanz <-conduct_ri(model_1 = held ~ 1, # restricted model
                            model_2 = held ~ genderint, # unrestricted model
                            assignment="genderint", declaration = dec,sharp_hypothesis = 0,data = iegendredTanz)
summary(anova_heldTanz) # p-value=.097
plot(anova_heldTanz)

###################################
####### Regression Analysis #######
###################################

# Independent Variables: 
# 1. gender of enumerator 
# 2. interaction term between gender of enumerator and gender of government official

# Dependent Variables: 
# 1. positive response
# 2. appointment held

# Covariates: 
# 1. status of the enumerator (local or foreign)
# 2. gender of the government official
# 3. social proof variable (control vs. social proof treatment)
# 4. modes of contact (email, unannounced visit, phone call)
# 5. status of government official-individual priority
# 6. interaction term between gender of enumerator and gender of government official
# 7. Country - fixed effect

# Cumulative Logistic Regression Analysis
fit_posresponse <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country + individualpriority + call + email + uvisit,
                         data = iegender,family = binomial) 
summary(fit_posresponse)

fit_posresponselm <- lm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country + individualpriority + call + email + uvisit,
                       data = iegender) 
summary(fit_posresponselm)
# why are these so different?

fit_held <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country + individualpriority + call + email + uvisit,
                       data = iegender,family = binomial) 
summary(fit_held)
# Interesting results in these two logistic regressions - being female increases the chances of getting a positive response and getting an appointment

# By country Logistic Regression Analysis
## India ##
fit_posresponseIndia <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + individualpriority + call + email + uvisit,
                       data = iegenderIndia,family = binomial) 
summary(fit_posresponseIndia)

fit_heldIndia <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority + call + email + uvisit,
                data = iegenderIndia,family = binomial) 
summary(fit_heldIndia)
# Gender not significant in India

## Tanzania ##
fit_posresponseTanz <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + individualpriority + call + email + uvisit,
                       data = iegenderTanz,family = binomial) 
summary(fit_posresponseTanz)

fit_heldTanz <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority + call + email + uvisit,
                data = iegenderTanz,family = binomial) 
summary(fit_heldTanz)
# gender significant and positive for appointments held, but not positive response

## Peru ##
fit_posresponsePeru <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority + call + email + uvisit,
                       data = iegenderPeru,family = binomial) 
summary(fit_posresponsePeru)

fit_heldPeru <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority + call + email + uvisit,
                data = iegenderPeru,family = binomial) 
summary(fit_heldPeru)
# Gender not significant in peru in this model

cor(iegender$call,as.numeric(iegender$Country))

#######################################################################################

# Test whether we should include methods of contact: call, email, uvisit

# test gender against methods of contact
ef <- glm(enumfemale~westerner + officialfemale + socialproof + email + call + uvisit + individualpriority, family=binomial,data=iegenderTanz)
summary(ef)

countmc <- lm(Country~email + call + uvisit,data=iegenderTanz)
summary(ef)

# test country against methods of contact
chisq.test(iegender$Country,iegender$call, correct=F)
chisq.test(iegender$Country,iegender$email, correct=F)
chisq.test(iegender$Country,iegender$uvisit, correct=F)

#######################################################################################

##########################################################
####### Regression Analysis without Contact Method #######
##########################################################

# Cumulative Logistic Regression Analysis
fit_posresponse <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country + individualpriority,
                       data = iegender,family = binomial) 
summary(fit_posresponse)

fit_posresponselm <- lm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country + individualpriority,
                        data = iegender) 
summary(fit_posresponselm)
# why are these so different?

fit_held <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country + individualpriority,
                data = iegender,family = binomial) 
summary(fit_held)
# Interesting results in these two logistic regressions - being female increases the chances of getting a positive response and getting an appointment

# By country Logistic Regression Analysis
## India ##
fit_posresponseIndia <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + individualpriority,
                            data = iegenderIndia,family = binomial) 
summary(fit_posresponseIndia)

fit_heldIndia <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority,
                     data = iegenderIndia,family = binomial) 
summary(fit_heldIndia)
# Gender not significant in India

## Tanzania ##
fit_posresponseTanz <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + individualpriority,
                           data = iegenderTanz,family = binomial) 
summary(fit_posresponseTanz)

fit_heldTanz <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority,
                    data = iegenderTanz,family = binomial) 
summary(fit_heldTanz)
# gender significant and positive for appointments held, but not positive response

## Peru ##
fit_posresponsePeru <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority,
                           data = iegenderPeru,family = binomial) 
summary(fit_posresponsePeru)

fit_heldPeru <- glm(held ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale  + individualpriority,
                    data = iegenderPeru,family = binomial) 
summary(fit_heldPeru)
# Gender not significant in peru in this model

##############################################
####### CACE - Appointments Held | Pos #######
##############################################

iegender.pos <- matrix(NA,14086,28)
for(i in 1:14086){
  for(j in 1:28){
  iegender.pos[i,j] <- ifelse(iegender$posresponse[i]=="1",iegender[i,j],NA)
  }
}

# Cumulative Difference of Means with Randomization Inference for Subset
iegender.pos <- subset(iegender, posresponse == 1)
# 779 data points with positive responses
with(iegender.pos, table(Country, enumfemale))
# Female Enumerators: 
# India - 255/398 (64.07%) 
# Peru - 124/238 (52.10%) 
# Tanzania - 64/143 (44.76%)
# Total - 321/779 (41.21%)

dec <- declare_ra(N = 779, m = 321)
ri_heldgpos <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                          assignment = "enumfemale", sharp_hypothesis = 0, data = iegender.pos)
summary(ri_heldgpos) # p-val=.419 - not significant
plot(ri_heldgpos)

# By Country Difference of Means with Randomization Inference for Subset
iegender.posIndia <- subset(iegenderIndia, posresponse == 1)
iegender.posTanz <- subset(iegenderTanz, posresponse == 1)
iegender.posPeru <- subset(iegenderPeru, posresponse == 1)

dec <- declare_ra(N = 398, m = 255)
ri_heldgposIndia <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                               assignment = "enumfemale", sharp_hypothesis = 0, data = iegender.posIndia)
summary(ri_heldgposIndia) # p-val=.404 - not significant
plot(ri_heldgposIndia)

dec <- declare_ra(N = 143, m = 64)
ri_heldgposTanz <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                              assignment = "enumfemale", sharp_hypothesis = 0, data = iegender.posTanz)
summary(ri_heldgposTanz) # p-val=.2 - not significant
plot(ri_heldgposTanz)

dec <- declare_ra(N = 238, m = 124)
ri_heldgposPeru <- conduct_ri(formula = held ~ enumfemale, declaration = dec, 
                          assignment = "enumfemale", sharp_hypothesis = 0, data = iegender.posPeru)
summary(ri_heldgposPeru) # p-val=.895 - not significant
plot(ri_heldgposPeru)

# None of these are significant

#############################################
####### Two-Stage Regression Analysis #######
#############################################

stageOne <- glm(posresponse ~ enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country,
                data = iegender,family = binomial)
xhat <- c(c(predict(stageOne)))
y <- rep(NA,length(xhat))
for(i in 1:length(xhat)) y[i] <- iegender$held[-is.na(iegender$enumfemale[i])]
stageTwo <- glm(y ~ xhat, family=binomial)

library(AER)
tsfit_posresponse <- ivreg(held ~ enumfemale + posresponse | enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + Country + individualpriority,
                           data = iegender)
summary(tsfit_posresponse) # positive response very significant. Gender is not.

tsfit_posresponseIndia <- ivreg(held ~ enumfemale + posresponse | enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + individualpriority,
                           data = iegenderIndia)
summary(tsfit_posresponseIndia) # positive response very significant. Gender is also significant - negative

tsfit_posresponseTanzania <- ivreg(held ~ enumfemale + posresponse | enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + individualpriority,
                           data = iegenderTanz)
summary(tsfit_posresponseTanzania) # positive response very significant. Gender is not.

tsfit_posresponsePeru <- ivreg(held ~ enumfemale + posresponse | enumfemale + westerner + officialfemale + socialproof + officialfemale*enumfemale + individualpriority,
                           data = iegenderPeru)
summary(tsfit_posresponsePeru) # positive response very significant. Gender is not.



### Logit regression: drop westerners and gender, and then run fixed effects for each enumerator
### Also include mode of contact (just to be true to our pre-analysis plan)

###########################################
####### Logit Regression with Names #######
###########################################

fit_posresponseNames <- glm(posresponse ~ Enumerator + officialfemale + socialproof + individualpriority,
                       data = iegender,family = binomial) 
summary(fit_posresponseNames)

fit_posresponseIndiaNames <- glm(posresponse ~ Enumerator + officialfemale + socialproof + individualpriority,
                            data = iegenderIndia,family = binomial) 
summary(fit_posresponseIndiaNames)

fit_posresponseTanzNames <- glm(posresponse ~ Enumerator + officialfemale + socialproof + individualpriority,
                           data = iegenderTanz,family = binomial) 
summary(fit_posresponseTanzNames)

fit_posresponsePeruNames <- glm(posresponse ~ Enumerator + officialfemale + socialproof  + individualpriority,
                           data = iegenderPeru,family = binomial) 
summary(fit_posresponsePeruNames)
