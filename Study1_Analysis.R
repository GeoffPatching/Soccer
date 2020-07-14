
#clear workspace
graphics.off()
rm(list=ls(all=TRUE))


###########################################################################################
#
#  TAKES ABOUT 2 HOURS TO RUN - and may disable your use of all other apps
#  You must have "Utilities.R" in your working directory after loading the data
#
###########################################################################################


# Select the directory with the data. Participant KB was removed, because KB consistently chose the left goal side over all the trials
setwd("~/Soccer/Study1_results_utanKB")


# Load the data from the directory
file_list = list.files()

index=1
for (file in file_list){
  
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=FALSE, quote="\"", comment.char="")
    index=index+1
    trials=c(1:nrow(temp_dataset))
    Subj= rep(index,nrow(temp_dataset))
    # get the response assignment "finger" used by participants.
    # For this to work response assignment "1" or "2" must be in the 6th column of the file name
    # To this end I renamed the following files...
    # 0126102014-May-2016.txt -> 0126102014-May-2016.txt, 
    # H120100102-May-2016.txt -> HA120100102-May-2016.txt
    # vb228100203-May-2016.txt -> vb228100203-May-2016.txt
    f=substr(file, start=6, stop=6)
    a=substr(file, start=4, stop=5)
    fing=rep(as.numeric(f),nrow(temp_dataset))
    age=rep(as.numeric(a),nrow(temp_dataset))
    temp_dataset=cbind(trials,Subj,age,fing,temp_dataset)
    temp_dataset=temp_dataset[-c(1:32),]
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=FALSE, quote="\"", comment.char="")
    Subj= rep(index,nrow(dataset))
    trials=c(1:nrow(dataset))
    # get the response assignment "finger" used by participants.
    # For this to work response assignment "1" or "2" must be in the 6th column of the file name
    # To this end I renamed the following files...
    # 0126102014-May-2016.txt -> 0126102014-May-2016.txt, 
    # H120100102-May-2016.txt -> HA120100102-May-2016.txt
    # vb228100203-May-2016.txt -> vb228100203-May-2016.txt
    f=substr(file, start=6, stop=6)
    a=substr(file, start=4, stop=5)
    fing=rep(as.numeric(f),nrow(dataset))
    age=rep(as.numeric(a),nrow(dataset))
    dataset=cbind(trials,Subj,age,fing,dataset)
    dataset=dataset[-c(1:32),]
  }
  
  
}

rm(list=setdiff(ls(), "dataset"))


# Select the working directory
setwd("~/Soccer")

# clean up the data and create transformed variables
# you must have Utilities.R in your working directory for this to work
source("Utilities.R")
data=tidy_up(dataset)

rm(dataset)

##################################################################
# Stan models
# Be patient, takes some hours to compute these models
# Note: may fail to compile if R session is not restarted
##################################################################

require(rethinking)

# leave 1 core available for other processes
nCores = parallel::detectCores() 
chains <- nCores <- nCores-1

# WARNING - HARDWARE DEPENDENT
# Increase the amount of memory available to R.
memory.limit(size = 7e+08)


# m1 regression of SRS on Logit P
# mean SRS and Logit P over all participants by condition

SummaryD=Summary(data)

m1 <- ulam(
  alist(
    logitP ~ normal( mu, sigma ),
    mu <- b0 + b1*SRS,
    c(b0,b1) ~ dnorm(0, 10),
    sigma ~ dunif(0 , 10)
  ) ,
  data=SummaryD, chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)

precis(m1,depth=2,prob=0.95)
post.m1 <- extract.samples(m1)



#m2 standard logistic regression of logit P on keeper and kicker
# m2 <- ulam(
#   alist(
#     Left ~ bernoulli( p ),
#     logit(p) <- b0 + b1*KeeperPer - b2*KickerPer,
#     c(b0,b1,b2) ~ dnorm(0,1)
#   ) ,
#   data=data, chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)
# 
# precis(m2,depth=2,prob=0.95,corr=TRUE)


#m3 logit P on Keeper, kicker , random intercepts and slopes
m3 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj] + b1[Subj]*KeeperPer - b2[Subj]*KickerPer,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m3,depth=3,prob=0.95)
post.m3 <- extract.samples(m3)


#m3.5 Logit P on Keeper, kicker + finger, random intercepts and slopes
# m3.5 <- ulam(
#   alist(
#     Left ~ bernoulli( p ),
#     logit(p) <- b0[Subj] + b1[Subj]*KeeperPer - b2[Subj]*KickerPer + b3*finger,
#     b0[Subj] ~ normal(b0_mu, b0_sigma),
#     b1[Subj] ~ normal(b1_mu, b1_sigma),
#     b2[Subj] ~ normal(b2_mu, b2_sigma),
#     b0_mu ~ normal(0,1),
#     b0_sigma ~ dcauchy(0,1),
#     b1_mu ~ normal(0,1),
#     b1_sigma ~ dcauchy(0,1),
#     b2_mu ~ normal(0,1),
#     b2_sigma ~ dcauchy(0,1),
#     b3 ~ normal(0,1)
#   ),
#   data=data,
#   constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0"),
#   warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)
# 
# precis(m3.5,depth=3,prob=0.95,corr=TRUE)

# compare(m3,m3.5)

#m4 SRS on Keeper, kicker, random intercepts and slopes
m4 <- ulam(
  alist(
    SRS ~ normal( mu, sigma ),
    mu <- b0[Subj] + b1[Subj]*KeeperPer - b2[Subj]*KickerPer,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1),
    sigma  ~ dcauchy(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0", sigma="lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m4,depth=3,prob=0.95)
post.m4 <- extract.samples(m4)


#m4 SRS on Keeper, kicker + finger, random intercepts and slopes

# m4.5 <- ulam(
#   alist(
#     SRS ~ normal( mu, sigma ),
#     mu <- b0[Subj] + b1[Subj]*KeeperPer - b2[Subj]*KickerPer + b3*finger,
#     b0[Subj] ~ normal(b0_mu, b0_sigma),
#     b1[Subj] ~ normal(b1_mu, b1_sigma),
#     b2[Subj] ~ normal(b2_mu, b2_sigma),
#     b0_mu ~ normal(0,1),
#     b0_sigma ~ dcauchy(0,1),
#     b1_mu ~ normal(0,1),
#     b1_sigma ~ dcauchy(0,1),
#     b2_mu ~ normal(0,1),
#     b2_sigma ~ dcauchy(0,1),
#     sigma  ~ dcauchy(0,1),
#     b3 ~ normal(0,1)
#   ),
#   data=data,
#   constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0", sigma="lower=0"),
#   warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)
# 
#  precis(m4.5,depth=3,prob=0.95,corr=TRUE)
# 
#  compare(m4,m4.5)


#m5 logit P on (Keeper-kicker)/2 and (Keeper+kicker)/2  , random intercepts and slopes
m5 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj] + b1[Subj]*KKPer_DiffAve + b2[Subj]*KKPer_SumAve,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m5,depth=3,prob=0.95,corr=TRUE)
post.m5 <- extract.samples(m5)



#m6 SRS  on (Keeper-kicker)/2 and (Keeper+kicker)/2  , random intercepts and slopes
m6 <- ulam(
  alist(
    SRS ~ normal( mu, sigma ),
    mu <- b0[Subj] + b1[Subj]*KKPer_DiffAve + b2[Subj]*KKPer_SumAve,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1),
    sigma  ~ dcauchy(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0", sigma="lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m6,depth=3,prob=0.95,corr=TRUE)
post.m6 <- extract.samples(m6)



####################################################
# Some fits by maximum likelihood
#####################################################

# SRS versus logit.P calculated over all participants
# This method avoids inf and 0
# windows()
# with(SummaryD, plot(logitP~SRS, xlab=expression(paste("Signed response speed (",italic(SRS),")")), ylab=expression(paste("Logit ", italic(P))),xlim=c(-1, 1), ylim=c(-1, 1),axes=T,bty="n"))
# abline(lm(logitP~SRS, data=SummaryD))
# text(0.5,-0.38, expression(paste("Adjusted ", italic(R)^2 == .97)))
# text(0.5,-0.5, expression(paste(hat(italic(y))," = 1.46",italic(x)," - 0.02")))
# summary(lm(logitP~SRS, data=SummaryD))
# 
# # Standared linear mixed effects via Maximum likelihood
# require(lme4)
# lm_m1=glmer(Left ~  KeeperPer + KickerPer + (1 + KeeperPer +  KickerPer | Subj), family=binomial(link = "logit"), data=data)
# summary(lm_m1)
# 
# lm_m1.5=glmer(Left ~  KeeperPer + KickerPer + finger + (1 + KeeperPer +  KickerPer | Subj), family=binomial(link = "logit"), data=data)
# summary(lm_m1.5)
# 
# anova(lm_m1,lm_m1.5)
# 
# lm_m2=lmer(SRS ~ ZkeeperPer + ZkickerPer + (1 + ZkeeperPer + ZkickerPer | Subj), data=data, REML=FALSE)
# summary(lm_m2)
# 
# lm_m2.5=lmer(SRS ~ ZkeeperPer + ZkickerPer + finger + (1 + ZkeeperPer + ZkickerPer | Subj), data=data, REML=FALSE)
# 
# summary(lm_m2.5)
# 
# anova(lm_m2,lm_m2.5)
# 
# lm_m3=glmer(Left ~  KKPer_DiffAve + KKPer_SumAve + (1 + KKPer_DiffAve +  KKPer_SumAve | Subj), family=binomial(link = "logit"), data=data)
# 
# summary(lm_m3)
# 
# lm_m4=lmer(SRS ~ KKPer_DiffAve + KKPer_SumAve + (1 + KKPer_DiffAve +  KKPer_SumAve | Subj), data=data, REML=FALSE)
# 
# summary(lm_m4)
# 
# 
# rm(lm_m1, lm_m1.5, lm_m2, lm_m2.5, lm_m3, lm_m4)
# 
# detach("package:lme4", unload=TRUE)


##################################################################
# Summary relations between SRS and Logit P
##################################################################

SRS_seq = seq(from=-1, to=1, by=0.1)
pred.data <- list(SRS=SRS_seq)
mu <- link(m1, data=pred.data)
mu.mean <- apply(mu,2,mean)
mu.HPDI <- apply(mu,2,HPDI, prob=0.95)


windows()
par(mar = c(5, 5, 2, 2))
plot(0, 0, type="n", xlab=expression(paste("Signed response speed (",italic(SRS),")")), ylab=expression(paste("Logit ", italic(P))),xlim=c(-1, 1), ylim=c(-1, 1),axes=T,bty="n", pch = 16)
lines(SRS_seq, mu.mean)
polygon(c(SRS_seq,rev(SRS_seq)),c(mu.HPDI[1,],rev(mu.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
with(SummaryD, points(SRS, logitP))


rm(mu.mean, SRS_seq, mu, mu.HPDI, pred.data)


###########################################
#counterfactual plots logit P and SRS on keeperPer, kicker
##########################################
windows(height=20,width=20)
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(2,2))

# Panel A conterfactual of logit P on keeperPer, kicker position held constant at 0
KeeperPer_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KeeperPer) post.m3$b0_mu + post.m3$b1_mu*KeeperPer
p <- sapply(KeeperPer_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="Goalkeeper position",ylab="Probability of left goal side selection", yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6","-5","-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(0, 1 , 0.1))
lines(KeeperPer_seq,p.mean)
polygon(c(KeeperPer_seq,rev(KeeperPer_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)

table=table(data$KeeperPer,data$Left)
PLeft=table[,2]/rowSums(table)
points(c(-4.9, -3.2, -1.6, 0, 1.6, 3.2, 4.9), PLeft)

mtext(text="(% left minus right difference in goalmouth area)",side=1,line=3.5)

# Panel B conterfactual of SRS on keeperPer, kicker position held constant at 0
KeeperPer_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KeeperPer) post.m4$b0_mu + post.m4$b1_mu*KeeperPer
p <- sapply(KeeperPer_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)

plot(0,0,type="n",xlab="Goalkeeper position",ylab=expression(paste("Signed response speed (",italic(SRS),")")), yaxt="n", ylim =c(-1,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6","-5","-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(-1, 1 , 0.2))
#mtext("Kicker aligned central to the veridical goalmouth", font=2, cex=1.5)
lines(KeeperPer_seq,p.mean)
polygon(c(KeeperPer_seq,rev(KeeperPer_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.0, lty=2)
segments(x0=0,y0=-1.2,x1=0,y1=1, lty=2)

SRS=with(data, tapply(SRS, KeeperPer, mean))
points(c(-4.9, -3.2, -1.6, 0, 1.6, 3.2, 4.9),SRS)

mtext(text="(% left minus right difference in goalmouth area)",side=1,line=3.5)

mtext("Kicker aligned central to the veridical goalmouth", outer = TRUE, line=-1, font=2, cex = 1.5)


# Panel c conterfactual logit P on kickerPer, keeper position held constant at 0
Kicker_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KickerPer) post.m3$b0_mu - post.m3$b2_mu*KickerPer
p <- sapply(Kicker_seq, p.link)
p.mean <- apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="Kicker position",ylab="Probability of left goal side selection", yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6","-5","-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(0, 1 , 0.1))
lines(Kicker_seq,p.mean)
polygon(c(Kicker_seq,rev(Kicker_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)

table=table(data$KickerPer,data$Left)
PLeft=table[,2]/rowSums(table)
points(c(-4.9, -3.2, -1.6, 0, 1.6, 3.2, 4.9), PLeft)

mtext(text="(% left minus right difference in goalmouth area)",side=1,line=3.5)


# Panel D conterfactual SRS on kickerPer, keeper position held constant at 0
Kicker_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KickerPer) post.m4$b0_mu - post.m4$b2_mu*KickerPer
p <- sapply(Kicker_seq, p.link)
p.mean <- apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)

plot(0,0,type="n",xlab="Kicker position",ylab=expression(paste("Signed response speed (",italic(SRS),")")), yaxt="n", ylim =c(-1,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6","-5","-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(-1, 1 , 0.2))
lines(Kicker_seq,p.mean)
polygon(c(Kicker_seq,rev(Kicker_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.0, lty=2)
segments(x0=0,y0=-1.2,x1=0,y1=1, lty=2)

SRS=with(data, tapply(SRS, KickerPer, mean))
points(c(-4.9, -3.2, -1.6, 0, 1.6, 3.2, 4.9),SRS)

mtext(text="(% left minus right difference in goalmouth area)",side=1,line=3.5)

mtext("Goalkeeper aligned central to the veridical goalmouth", outer = TRUE, font=2, line=-29, cex=1.5)

###########################################
#counterfactual plots logit P and SRS on (keeperPer-kicker)/2 and (Keeper+kicker)/2
##########################################
windows(height=20,width=20)
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(2,2))

# Panel A conterfactual logit P on (keeperPer-kicker)/2, (Keeper+kicker)/2 position held constant at 0
KK_seq <- seq(from=-4, to=4, by=0.5)
p.link <- function(KKPer_DiffAve) post.m5$b0_mu + post.m5$b1_mu*KKPer_DiffAve
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="(Goalkeeper position - Kicker position) / 2",ylab="Probability of left goal side selection", yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-4,4),bty="n")
axis(1, at=-4:4, labels=c("-4","-3","-2","-1","0","+1","+2","+3","+4"))
axis(2, at=seq(0, 1 , 0.1))
#mtext("Goalkeeper position + Kicker position = 0", font=2, cex=1.5)
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)

table=table(data$KKPer_DiffAve,data$Left)
PLeft=table[,2]/rowSums(table)
points(c(-2.4, -0.8, 0.8, 2.4), PLeft)

#mtext(text="% left goalmouth area defined by the two players relative difference from central",side=1,line=3.5)


# Panel B conterfactual SRS on (keeperPer-kicker)/2, (Keeper+kicker)/2 position held constant at 0
KK_seq <- seq(from=-4, to=4, by=0.5)
p.link <- function(KKPer_DiffAve) post.m6$b0_mu + post.m6$b1_mu*KKPer_DiffAve
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)


plot(0,0,type="n",xlab="(Goalkeeper position - Kicker position) / 2",ylab=expression(paste("Signed response speed (",italic(SRS),")")), yaxt="n", ylim =c(-1,1), xaxt="n", xlim=c(-4,4),bty="n")
axis(1, at=-4:4, labels=c("-4","-3","-2","-1","0","+1","+2","+3","+4"))
axis(2, at=seq(-1, 1 , 0.2))
#mtext("Goalkeeper position + Kicker position = 0", font=2, cex=1.5)
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.0, lty=2)
segments(x0=0,y0=-1.2,x1=0,y1=1, lty=2)

SRS=with(data, tapply(SRS, KKPer_DiffAve, mean))
points(c(-2.4, -0.8, 0.8, 2.4),SRS)


#mtext(text="% left goalmouth area defined by the two players relative difference from central",side=1,line=3.5)

mtext("Goalkeeper position + Kicker position = 0", outer = TRUE, line=-1, font=2, cex=1.5)


# Panel c conterfactual logit P on (Keeper+kicker)/2, (keeperPer-kicker)/2 position held constant at 0
KK_seq <- seq(from=-4, to=4, by=0.5)
p.link <- function(KKPer_SumAve) post.m5$b0_mu + post.m5$b2_mu*KKPer_SumAve
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="(Goalkeeper position + Kicker position) / 2",ylab="Probability of left goal side selection", yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-4,4),bty="n")
axis(1, at=-4:4, labels=c("-4","-3","-2","-1","0","+1","+2","+3","+4"))
axis(2, at=seq(0, 1 , 0.1))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)

table=table(data$KKPer_SumAve,data$Left)
PLeft=table[,2]/rowSums(table)
points(c(-2.4, -0.8, 0.8, 2.4), PLeft)


#mtext(text="% left goalmouth area defined by the two players average position from central",side=1,line=3.5)


# Panel D conterfactual SRS on (Keeper+kicker)/2, (keeperPer-kicker)/2 position held constant at 0
KK_seq <- seq(from=-4, to=4, by=0.5)
p.link <- function(KKPer_SumAve) post.m6$b0_mu + post.m6$b2_mu*KKPer_SumAve
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)


plot(0,0,type="n",xlab="(Goalkeeper position + Kicker position) / 2",ylab=expression(paste("Signed response speed (",italic(SRS),")")), yaxt="n", ylim =c(-1,1), xaxt="n", xlim=c(-4,4),bty="n")
axis(1, at=-4:4, labels=c("-4","-3","-2","-1","0","+1","+2","+3","+4"))
axis(2, at=seq(-1, 1 , 0.2))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.0, lty=2)
segments(x0=0,y0=-1.2,x1=0,y1=1, lty=2)

SRS=with(data, tapply(SRS, KKPer_SumAve, mean))
points(c(-2.4, -0.8, 0.8, 2.4),SRS)

#mtext(text="% left goalmouth area defined by the two players average position from central",side=1,line=3.5)

mtext("Goalkeeper position - Kicker position = 0", outer = TRUE, font=2, line=-29, cex=1.5)

#######################################################################
# That's it
# Have fun modifying the script to try out different possibilities.
#######################################################################



