# Data Loading
bank <- read.csv("/Users/Chenyu/Downloads/bank-additional-full.csv", header=TRUE, stringsAsFactors=FALSE)

bank$job <- as.factor(bank$job)
bank$marital <- as.factor(bank$marital)
bank$education <- as.factor(bank$education)
bank$housing <- as.factor(bank$housing)
bank$loan <- as.factor(bank$loan)
bank$contact <- as.factor(bank$contact)
bank$month <- as.factor(bank$month)
bank$day_of_week <- as.factor(bank$day_of_week)
bank$poutcome <- as.factor(bank$poutcome)
bank$y <- as.factor(bank$y)
bank$ynum <- as.numeric(bank$y) - 1

# install.packages("naniar")
library(naniar)
bankna <- replace_with_na(bank, replace = list(job ="unknown", marital = "unknown", education = "unknown",
                                               housing = "unknown", loan = "unknown", poutcome = "unknown"))
na_count <-sapply(bankna, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
bank = subset(bankna, select = -c(default,duration))
summary(bank)

library(mice)
# Visualize Missing Data
md.pattern(bank)
# install.packages("VIM")
library(VIM);library(lattice)
aggr(bank,col=c("lightblue3","darkred"),numbers=TRUE,sortVars=TRUE,labels=names(bank),cex.axis=.7,gap=3,ylab=c("Proportion missing","Missingness pattern"))
marginplot(bank[,c("diameter.","age")],col=c("lightblue3","darkred"),cex.numbers = 1.2,pch=19)
# Imputation
bank_imp <- mice(bank,m=5,print=F)

library(ggplot2)
d1 <- complete(bank_imp,1);d1
apply(table(d1[,c("y","education")])/sum(table(d1[,c("y","education")])),2,function(x) x/sum(x))
apply(table(d1[,c("y","housing")])/sum(table(d1[,c("y","housing")])),2,function(x) x/sum(x))
apply(table(d1[,c("y","loan")])/sum(table(d1[,c("y","loan")])),2,function(x) x/sum(x))
d1$oldconsumer[d1$pdays == 999] <- 0
d1$oldconsumer[d1$pdays < 999] <- 1
table(d1$oldconsumer)
d1$newage[d1$age <=30] <- 0
d1$newage[d1$age >30 & d1$age < 60] <- 1
d1$newage[d1$age >=60] <- 2
d1$newage <- as.factor(d1$newage)
table(d1$newage)

d2 <- complete(bank_imp,2);d2
apply(table(d2[,c("y","education")])/sum(table(d2[,c("y","education")])),2,function(x) x/sum(x))
apply(table(d2[,c("y","housing")])/sum(table(d2[,c("y","housing")])),2,function(x) x/sum(x))
apply(table(d2[,c("y","loan")])/sum(table(d2[,c("y","loan")])),2,function(x) x/sum(x))

bank_obs <- na.omit(bank)
apply(table(bank_obs[,c("y","education")])/sum(table(bank_obs[,c("y","education")])),2,function(x) x/sum(x))
apply(table(bank_obs[,c("y","housing")])/sum(table(bank_obs[,c("y","housing")])),2,function(x) x/sum(x))
apply(table(bank_obs[,c("y","loan")])/sum(table(bank_obs[,c("y","loan")])),2,function(x) x/sum(x))

library(MatchIt) #for propensity score matching
library(cobalt)
library(tidyverse)
library(knitr)
library(GGally)
library(xtable)
library(arm)
library(pROC)
library(e1071)
library(caret)
library(rms)
require(gridExtra)

# EDA

# age (none / people are less likely to purchase the product during their mid-ages (30s-50s), compared to their 20s and 60s)
table(bank$age)
ggplot(bank,aes(x=y, y=age, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Age vs Termed Deposit",
       x="Termed Deposit?",y="Age") + 
  theme_classic() + theme(legend.position="none")
binnedplot(y=bank$ynum,bank$age,xlab="Age",ylim=c(0,1),col.pts="navy",
           ylab ="Termed Deposit?",main="Binned Age and Term Deposits",
           col.int="white")


# job (retired - highest)
table(bank$job)
tjob <- apply(table(bank[,c("y","job")])/sum(table(bank[,c("y","job")])),2,function(x) x/sum(x))
xtable(tjob)
plot(0:11,tapply(bank$ynum, bank$job, mean),col='blue4',pch=10)

# marital (similar)
table(bank$marital)
tmarital <- apply(table(bank[,c("y","marital")])/sum(table(bank[,c("y","marital")])),2,function(x) x/sum(x))
xtable(tmarital)
plot(0:3,tapply(bank$ynum, bank$marital, mean),col='blue4',pch=10)

# education (high school - highest)
table(bank$education)
teducation <- apply(table(bank[,c("y","education")])/sum(table(bank[,c("y","education")])),2,function(x) x/sum(x))
xtable(teducation)
plot(0:7,tapply(bank$ynum, bank$education, mean),col='blue4',pch=10)

# housing (similar)
table(bank$housing)
thousing <- apply(table(bank[,c("y","housing")])/sum(table(bank[,c("y","housing")])),2,function(x) x/sum(x))
xtable(thousing)
plot(0:2,tapply(bank$ynum, bank$housing, mean),col='blue4',pch=10)

# loan (similar)
table(bank$loan)
tloan <- apply(table(bank[,c("y","loan")])/sum(table(bank[,c("y","loan")])),2,function(x) x/sum(x))
xtable(tloan)
plot(0:2,tapply(bank$ynum, bank$loan, mean),col='blue4',pch=10)

# contact (cellular slightly higher)
tcontact <- apply(table(bank[,c("y","contact")])/sum(table(bank[,c("y","contact")])),2,function(x) x/sum(x))
xtable(tcontact)
plot(0:1,tapply(bank$ynum, bank$contact, mean),col='blue4',pch=10)

# month (apr, aug, jul, jun, may, nov)
tmonth <- apply(table(bank[,c("y","month")])/sum(table(bank[,c("y","month")])),2,function(x) x/sum(x))
xtable(tmonth)
plot(0:9,tapply(bank$ynum, bank$month, mean),col='blue4',pch=10)

# day (similar)
tday <- apply(table(bank[,c("y","day_of_week")])/sum(table(bank[,c("y","day_of_week")])),2,function(x) x/sum(x))
xtable(tday)
plot(0:4,tapply(bank$ynum, bank$day_of_week, mean),col='blue4',pch=10)

# duration (positive) ? don't include
ggplot(bank,aes(x=y, y=duration, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Duration vs Termed Deposit",
       x="Termed Deposit?",y="Duration") + 
  theme_classic() + theme(legend.position="none")
binnedplot(y=bank$ynum,bank$duration,xlab="Duration",ylim=c(0,1),col.pts="navy",
           ylab ="Termed Deposit?",main="Binned Duration and Termed Deposit",
           col.int="white")
cor(bank$ynum,bank$duration)
regduration <- glm(y~duration, data=bank, family = binomial)
summary(regduration)

# campaign (none / decreasing over the number of contacts)
ggplot(bank,aes(x=y, y=campaign, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Campaign vs Termed Deposit",
       x="Termed Deposit?",y="Campaign") + 
  theme_classic() + theme(legend.position="none")
binnedplot(y=bank$ynum,bank$campaign,xlab="Campaign",ylim=c(0,1),col.pts="navy",
           ylab ="Term Deposits?",main="Binned Number of Contacts During this Campaign and Term Deposits",
           col.int="white")

# previous
ggplot(bank,aes(x=y, y=previous, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Contacts for this Compaign vs Termed Deposit",
       x="Termed Deposit?",y="Contacts for this Compaign") + 
  theme_classic() + theme(legend.position="none")
bank$yescontact[bank$previous == 0] <- 0
bank$yescontact[bank$previous != 0] <- 1
tyescontact <- apply(table(bank[,c("y","yescontact")])/sum(table(bank[,c("y","yescontact")])),2,function(x) x/sum(x))
xtable(tyescontact)
binnedplot(y=bank$ynum,bank$previous,xlab="Contacts before this Campaign",ylim=c(0,1),col.pts="navy",
           ylab ="Term Deposit?",main="Binned Contacts before this Campaign and Term Deposits",
           col.int="white")
# plot(0:1,tapply(bank$ynum, bank$yescontact, mean),col='blue4',pch=10)

# poutcome: consumers who have previously bought termed deposits are more likely to purchase again
table(bank$poutcome)
tpoutcome <- apply(table(bank[,c("y","poutcome")])/sum(table(bank[,c("y","poutcome")])),2,function(x) x/sum(x))
xtable(tpoutcome)
plot(0:2,tapply(bank$ynum, bank$poutcome, mean),col='blue4',pch=10)

# employment variation rate (negative) [drop - high correlation with CPI and Euribor]
ggplot(bank,aes(x=y, y=emp.var.rate, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Employment Variation Rate vs Termed Deposit",
       x="Termed Deposit?",y="Employment Variation Rate") + 
  theme_classic() + theme(legend.position="none")
binnedplot(y=bank$ynum,bank$emp.var.rate,xlab="Employment Variation Rate",ylim=c(0,1),col.pts="navy",
           ylab ="Term Deposit?",main="Binned Employment Variation Rate and Term Deposits",
           col.int="white")

# consumer price index (negative)
ggplot(bank,aes(x=y, y=cons.price.idx, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Consumer Price Index vs Term Deposits",
       x="Term Deposits?",y="Consumer Price Index") + 
  theme_classic() + theme(legend.position="none")
binnedplot(y=bank$ynum,bank$cons.price.idx,xlab="Consumer Price Index",ylim=c(0,1),col.pts="navy",
           ylab ="Term Deposits?",main="Binned Consumer Price Index and Term Deposits",
           col.int="white")

# CCI (positive)
ggplot(bank,aes(x=y, y=cons.conf.idx, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="CCI vs Term Deposits",
       x="Term Deposits?",y="CCI") + 
  theme_classic() + theme(legend.position="none")
binnedplot(y=bank$ynum,bank$cons.conf.idx,xlab="CCI",ylim=c(0,1),col.pts="navy",
           ylab ="Term Deposits?",main="Binned CCI and Term Deposits",
           col.int="white")

# Euribor (negative)
ggplot(bank,aes(x=y, y=euribor3m, fill=y)) +
  geom_boxplot() + coord_flip() +
  scale_fill_brewer(palette="Reds") +
  labs(title="Euribor 3 Month Rates vs Term Deposits",
       x="Term Deposits?",y="Euribor 3 Month Rates") + 
  theme_classic() + theme(legend.position="none")
binnedplot(y=bank$ynum,bank$euribor3m,xlab="Euribor 3 Month Rates",ylim=c(0,1),col.pts="navy",
           ylab ="Term Deposits?",main="Binned Euribor 3 Month Rates and Term Deposits",
           col.int="white")
summary(d1$euribor3m)
d1$loweuri[d1$euribor3m<=3]<-1
d1$loweuri[d1$euribor3m>3]<-0
d1$loweuri<-as.factor(d1$loweuri)

# correlation: employment variation rate, consumer price index, CCI, Euribor
cor(bank$emp.var.rate, bank$cons.conf.idx)
cor(bank$cons.price.idx, bank$emp.var.rate)
cor(bank$emp.var.rate, bank$euribor3m)
cor(bank$cons.price.idx, bank$cons.conf.idx)
cor(bank$cons.price.idx, bank$euribor3m)
cor(bank$cons.conf.idx, bank$euribor3m)

# interaction
# age * poutcome
binnedplot(d1$agec[d1$poutcome=="failure"], 
           y=d1$ynum[d1$poutcome=="failure"], 
           xlab = "Age-centered", ylab = "Termed Deposit", main = "Binned Age-centered and Termed Deposit (failure)") 
binnedplot(d1$agec[d1$poutcome=="nonexistent"], 
           y=d1$ynum[d1$poutcome=="nonexistent"], 
           xlab = "Age-centered", ylab = "Termed Deposit", main = "Binned Age-centered and Termed Deposit (nonexist)") 
binnedplot(d1$agec[d1$poutcome=="success"], 
           y=d1$ynum[d1$poutcome=="success"], 
           xlab = "Age-centered", ylab = "Termed Deposit", main = "Binned Age-centered and Termed Deposit (success)") 

# newage:poutcome
plot(0:2,tapply(d1$ynum[d1$poutcome=="failure"], d1$newage[d1$poutcome=="failure"], mean),col='blue4',pch=10,ann=FALSE,ylim=c(0,1),cex=1.5)
par(new=TRUE)
plot(0:2,tapply(d1$ynum[d1$poutcome=="nonexistent"], d1$newage[d1$poutcome=="nonexistent"], mean),col='red4',pch=10,ann=FALSE,ylim=c(0,1),cex=1.5)
par(new=TRUE)
plot(0:2,tapply(d1$ynum[d1$poutcome=="success"], d1$newage[d1$poutcome=="success"], mean),col='black',pch=10,ann=FALSE,ylim=c(0,1),cex=1.5)
title(xlab = "Age",ylab = 'Prob of Purchasing a Term Deposit')
legend("topright",pch=c(10,10),legend=c("failure","nonexistent","success"),col=c("blue4","red4","black"),bty="n")

# poutcome:confc
binnedplot(d1$confc[d1$poutcome=="failure"], 
           y=d1$ynum[d1$poutcome=="failure"], 
           xlab = "CCI-centered", ylab = "Term Deposits", main = "Binned CCI-centered and Term Deposits (failure)") 
binnedplot(d1$confc[d1$poutcome=="nonexistent"], 
           y=d1$ynum[d1$poutcome=="nonexistent"], 
           xlab = "CCI-centered", ylab = "Term Deposits", main = "Binned CCI-centered and Term Deposits (nonexist)") 
binnedplot(d1$confc[d1$poutcome=="success"], 
           y=d1$ynum[d1$poutcome=="success"], 
           xlab = "CCI-centered", ylab = "Term Deposits", main = "Binned CCI-centered and Term Deposits (success)") 

# poutcome:loweuri
plot(0:1,tapply(d1$ynum[d1$poutcome=="failure"], d1$loweuri[d1$poutcome=="failure"], mean),col='blue4',pch=10,ann=FALSE,ylim=c(0,1),cex=1.5)
par(new=TRUE)
plot(0:1,tapply(d1$ynum[d1$poutcome=="nonexistent"], d1$loweuri[d1$poutcome=="nonexistent"], mean),col='red4',pch=10,ann=FALSE,ylim=c(0,1),cex=1.5)
par(new=TRUE)
plot(0:1,tapply(d1$ynum[d1$poutcome=="success"], d1$loweuri[d1$poutcome=="success"], mean),col='black',pch=10,ann=FALSE,ylim=c(0,1),cex=1.5)
title(xlab = "Euribor",ylab = 'Prob of Purchasing a Term Deposit')
legend("topright",pch=c(10,10),legend=c("failure","nonexistent","success"),col=c("blue4","red4","black"),bty="n")


# centering
bank$agec <- bank$age-mean(bank$age)
bank$varc <- bank$emp.var.rate + mean(bank$emp.var.rate)
bank$pricec <- bank$cons.price.idx - mean(bank$cons.price.idx)
bank$confc <- bank$cons.conf.idx - mean(bank$cons.conf.idx)
bank$euric <- bank$euribor3m - mean(bank$euribor3m)

d1$agec <- d1$age-mean(d1$age)
d1$varc <- d1$emp.var.rate + mean(d1$emp.var.rate)
d1$pricec <- d1$cons.price.idx - mean(d1$cons.price.idx)
d1$confc <- d1$cons.conf.idx - mean(d1$cons.conf.idx)d
d1$euric <- d1$euribor3m - mean(d1$euribor3m)

# model 1
model1 <- glm(ynum~ newage + job + marital + education + housing + loan + contact
              + month + day_of_week + campaign + previous + poutcome
              + varc + pricec + confc + loweuri, data = d1, family = binomial)
summary(model1)

rawresid1 <- residuals(model1,"resp")
binnedplot(x=fitted(model1),y=rawresid1,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$campaign,y=rawresid1,xlab="Number of Contacts During this Campaign",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$previous,y=rawresid1,xlab="Contacts for this Campaign",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$varc,y=rawresid1,xlab="employment variation rate centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$pricec,y=rawresid1,xlab="CPI centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$confc,y=rawresid1,xlab="CCI centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")

# Model validation
# multicollinearity
vif(model1)
#let's do the confusion matrix with .5 threshold
Conf_mat <- confusionMatrix(as.factor(ifelse(fitted(model1) >=  mean(d1$ynum), "1","0")),
                            as.factor(d1$ynum),positive = "1")
Conf_mat$table
Conf_mat$overall["Accuracy"];
Conf_mat$byClass[c("Sensitivity","Specificity")] #True positive rate and True negative rate
#Maybe we can try to increase that accuracy.
roc(bank$ynum,fitted(model1),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3")

model2 <- glm(ynum~ newage + job + marital + education + housing + loan + previous
              + contact + month + day_of_week + campaign + poutcome + confc + loweuri,
              data = d1, family = binomial)
summary(model2)
vif(model2)

# interaction
model3 <- glm(ynum~ newage + job + marital + education + housing + loan + previous
              + contact + month + day_of_week + campaign + poutcome
              + confc + pricec + loweuri + poutcome:confc, data = d1, family = binomial)
summary(model3)
vif(model3)

# Stepwise
model0 <- glm(ynum~ 1, data = d1, family = binomial)
model_stepwiseaic <- step(model0,scope=formula(model3),direction="both", trace=0)
model_stepwiseaic$call
Model_aic <- glm(ynum~ loweuri + month + poutcome + confc + pricec + 
                   contact + day_of_week + job + campaign + newage + education + 
                   poutcome:confc,
                 data = d1, family = binomial)
summary(Model_aic)
anova(model3, Model_aic, test= "Chisq")
# no difference: pick Model_aic

vif(Model_aic)
rawresid2 <- residuals(Model_aic,"resp")
binnedplot(x=fitted(Model_aic),y=rawresid2,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$previous,y=rawresid2,xlab="Contacts for this Campaign",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$confc,y=rawresid2,xlab="CCI centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$pricec,y=rawresid2,xlab="CPI centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$campaign,y=rawresid2,xlab="Number of Contacts During this Campaign",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")

Conf_mat <- confusionMatrix(as.factor(ifelse(fitted(Model_aic) >=  mean(d1$ynum), "1","0")),
                            as.factor(d1$ynum),positive = "1")
Conf_mat$table
Conf_mat$overall["Accuracy"];
Conf_mat$byClass[c("Sensitivity","Specificity")] #True positive rate and True negative rate
#Maybe we can try to increase that accuracy.
roc(d1$ynum,fitted(Model_aic),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3")

# CPI
model4 <- glm(ynum~ poutcome + month + contact + newage + pricec + 
                job + campaign + day_of_week + education + marital+poutcome:pricec,
              data = d1, family = binomial)
vif(model4)
rawresid4 <- residuals(model4,"resp")
binnedplot(x=fitted(model4),y=rawresid4,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$pricec,y=rawresid4,xlab="CPI centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")

model5 <- glm(ynum~ poutcome + month + contact + newage + confc + pricec + loweuri+ previous +
                job + campaign + day_of_week + education + marital+ poutcome:confc,
              data = d1, family = binomial)
summary(model5)
vif(model5)
rawresid5 <- residuals(model5,"resp")
binnedplot(x=fitted(model5),y=rawresid5,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$confc,y=rawresid5,xlab="CCI centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$pricec,y=rawresid5,xlab="CPI centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$campaign,y=rawresid5,xlab="Number of Contacts During this Campaign",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")


# Hierachical Model
multimodel5 <- glmer(formula = ynum~poutcome + (1|month) + contact + newage + confc + pricec + loweuri + previous +
                       job + campaign + day_of_week + education + marital + poutcome:confc, family = binomial(link="logit"), 
                     data = d1)
summary(multimodel5)
library(sjPlot)
tab_model(multimodel5)
dotplot(ranef(multimodel5, condVar=TRUE))

rawresid5 <- residuals(multimodel5,"resp")
binnedplot(x=fitted(multimodel5),y=rawresid5,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$previous,y=rawresid5,xlab="Number of Contacts Before this Campaign",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$confc,y=rawresid5,xlab="CCI-centered",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
binnedplot(x=d1$campaign,y=rawresid5,xlab="Number of Contacts During this Campaign",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
invisible(roc(d1$ynum,fitted(multimodel5),plot=T,print.thres="best",legacy.axes=T,print.auc =T,col="red3"))


(ranef(multimodel5)$month)["apr",]

Conf_mat <- confusionMatrix(as.factor(ifelse(fitted(multimodel5) >=  mean(d1$ynum), "1","0")),
                            as.factor(d1$ynum),positive = "1")
Conf_mat$table
Conf_mat$overall["Accuracy"];
Conf_mat$byClass[c("Sensitivity","Specificity")] #True positive rate and True negative rate
#Maybe we can try to increase that accuracy.
roc(d1$ynum,fitted(multimodel5),plot=T,print.thres="best",legacy.axes=T,
    print.auc =T,col="red3")

# Eliminate Outlier: May & March
newd1 <- d1[!(d1$month=="may"),]
newd1 <- newd1[!(newd1$month=="mar"),]
newmodel <- glmer(formula = ynum~poutcome + (1|month) + contact + newage + confc + pricec + loweuri + previous +
                    job + campaign + day_of_week + education + marital + poutcome:confc, family = binomial(link="logit"), 
                  data = newd1)
summary(newmodel)
newrawresid <- residuals(newmodel,"resp")
binnedplot(x=fitted(newmodel),y=newrawresid,xlab="Pred. probabilities",
           col.int="red4",ylab="Avg. residuals",main="Binned residual plot",col.pts="navy")
