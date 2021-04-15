
#clear workspace
graphics.off()
rm(list=ls(all=TRUE))


###########################################################################################
#
#  Takes some hours to run
#  You may save time by running fewer warmup trials / chains 
#  You must have "Study1_DataPrep.R" in your working directory after loading the data
#
###########################################################################################


# Select the directory with the data. Participant KB was removed, because KB consistently chose the left goal side over all the trials
setwd("~/Soccer/Study1_Data_utanKB")


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
    s=substr(file, start=3, stop=3)
    fing=rep(as.numeric(f),nrow(temp_dataset))
    age=rep(as.numeric(a),nrow(temp_dataset))
    sex=rep(as.numeric(s),nrow(temp_dataset))
    temp_dataset=cbind(trials,Subj,age,fing,sex,temp_dataset)
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
    s=substr(file, start=3, stop=3)
    fing=rep(as.numeric(f),nrow(dataset))
    age=rep(as.numeric(a),nrow(dataset))
    sex=rep(as.numeric(s),nrow(dataset))
    dataset=cbind(trials,Subj,age,fing,sex,dataset)
    dataset=dataset[-c(1:32),]
  }
  
  
}

rm(list=setdiff(ls(), "dataset"))


# clean up the data and create transformed variables
# you must have Study1_DataPrep.R in your working directory for this to work

# Select the working directory
setwd("~/Soccer")
source("Study1_DataPrep.R")
data=tidy_up(dataset)

rm(dataset,Summary,SummaryBySubj,tidy_up)
detach(package:psych)

load("models.RData")

##################################################################
# Stan models
# Be patient, takes about 30 mins to compute these models
##################################################################

require(rethinking)

# leave 1 core available for other processes
nCores = parallel::detectCores() 
chains <- nCores <- nCores-1


#m1 logit P on Keeper, kicker , random intercepts and slopes
m1 <- ulam(
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

precis(m1,depth=3,prob=0.95)
post.m1 <- extract.samples(m1)


###########################################
#counterfactual plots logit P on keeper, kicker
##########################################
windows(height=15,width=30) #use X11 for mac
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))

# Panel A conterfactual of logit P on keeper, kicker position held constant at 0
KeeperPer_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KeeperPer) post.m1$b0_mu + post.m1$b1_mu*KeeperPer
p <- sapply(KeeperPer_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="Goalkeeper position",ylab=expression(paste(italic(P),"(left goal side)")), yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-6,6),bty="n")
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
mtext("Kicker aligned central to the veridical goalmouth", side=3, line=1, font=2, cex = 1.5)


# Panel B conterfactual logit P on kicker, keeper position held constant at 0
Kicker_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KickerPer) post.m1$b0_mu - post.m1$b2_mu*KickerPer
p <- sapply(Kicker_seq, p.link)
p.mean <- apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="Kicker position",ylab=expression(paste(italic(P),"(left goal side)")), yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-6,6),bty="n")
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
mtext("Goalkeeper aligned central to the veridical goalmouth", side=3, font=2, line=1, cex=1.5)

############################################################################
############################################################################

#m1a logit P on (Keeper-kicker) and (Keeper+kicker), random intercepts and slopes
m1a <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj] + b1[Subj]*KKPerDiff + b2[Subj]*KKPerSum,
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

precis(m1a,depth=3,prob=0.95,corr=TRUE)
post.m1a <- extract.samples(m1a)


###########################################
#counterfactual plots logit P and SRS on (keeper-kicker) and (Keeper+kicker)
##########################################
windows(height=15,width=30) #use X11 for mac
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))

# Panel A conterfactual logit P on (keeper-kicker), (Keeper+kicker) position held constant at 0
KK_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KKPerDiff) post.m1a$b0_mu + post.m1a$b1_mu*KKPerDiff
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="Goalkeeper position - Kicker position",ylab=expression(paste(italic(P),"(left goal side)")), yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6", "-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(0, 1 , 0.1))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)

table=table(data$KKPerDiff,data$Left)
PLeft=table[,2]/rowSums(table)
points(c(-4.8, -1.6, 1.6, 4.8), PLeft)
mtext("Goalkeeper position + Kicker position = 0", side=3, line=1, font=2, cex=1.5)

# Panel B conterfactual logit P on (Keeper+kicker), (keeper-kicker) position held constant at 0
KK_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KKPerSum) post.m1a$b0_mu + post.m1a$b2_mu*KKPerSum
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",xlab="Goalkeeper position + Kicker position",ylab=expression(paste(italic(P),"(left goal side)")), yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6", "-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(0, 1 , 0.1))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)

table=table(data$KKPerSum,data$Left)
PLeft=table[,2]/rowSums(table)
points(c(-4.8, -1.6, 1.6, 4.8), PLeft)
mtext("Goalkeeper position - Kicker position = 0", side=3, line=1, font=2, cex=1.5)



#######################################################################
# Other analyses of interest and model comparisons
#######################################################################
require(rethinking)

# leave 1 core available for other processes
nCores = parallel::detectCores() 
chains <- nCores <- nCores-1

#m2 logit P on finger (response assignment)
m2 <- ulam(
  alist(
    Left ~ binomial(1, p ),     
    logit(p) <- a[finger],
    a[finger] ~ normal(0, 1)          
  ),
  data=data, warmup=500 , iter=1000 , chains=1 , cores=nCores)

precis(m2,depth=3,prob=0.95)
post <- extract.samples(m2)
post$diff_fm <- post$b0 - post$b1
precis( post , depth=3,hist=FALSE,prob=0.95 )
rm(m2,post)

#m3 logit P on sex (male = 1, female = 2)
m3 <- ulam(
  alist(
    Left ~ binomial(1, p ),     
    logit(p) <- a[sex],
    a[sex] ~ normal(0, 1)          
  ),
  data=data, warmup=500 , iter=1000 , chains=1 , cores=nCores)

precis(m3,depth=3,prob=0.95)
post <- extract.samples(m3)
post$diff_fm <- post$b0 - post$b1
precis( post , depth=3,hist=FALSE,prob=0.95 )

rm(m3,post)

#m4 logit P intercept only
m4 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj],
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m4,depth=3,prob=0.95)

#m5 logit P intercept + keeperpos
m5 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj] + b1[Subj]*KeeperPer,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m5,depth=3,prob=0.95)

compare(m4,m5) #intercept only vs intercept + keeper
compare(m5,m1) #intercept + keeper vs intercept + keeper + kicker


#m6 logit P intercept + (Keeper-kicker)
m6 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj] + b1[Subj]*KKPerDiff,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m6,depth=3,prob=0.95)

compare(m4,m6) #intercept only vs intercept + keeperdiff
compare(m6,m1a) #intercept + keeperdiff vs intercept + keeperdiff +keepersum



#m7 Logit P on Keeper, kicker + finger, random intercepts and slopes
m7 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- a[finger]+b0[Subj] + b1[Subj]*KeeperPer - b2[Subj]*KickerPer,
    a[finger] ~ normal(0, 1),
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1),
    b3 ~ normal(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m7,depth=3,prob=0.95,corr=TRUE)


compare(m1,m7) #intercept + keeper + kicker vs. intercept + keeper + kicker + finger


#m8 Logit P on Keeper, kicker + sex, random intercepts and slopes
m8 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- a[sex]+b0[Subj] + b1[Subj]*KeeperPer - b2[Subj]*KickerPer,
    a[sex] ~ normal(0, 1),
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1),
    b3 ~ normal(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m8,depth=3,prob=0.95,corr=TRUE)

compare(m1,m8) #intercept + keeper + kicker vs. intercept + keeper + kicker + sex

m9 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj] + b1[Subj]*KKPerDiff + b2[Subj]*KKPerSum + b3*finger,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1),
    b3 ~ normal(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)


compare(m1a,m9)


m10 <- ulam(
  alist(
    Left ~ bernoulli( p ),
    logit(p) <- b0[Subj] + b1[Subj]*KKPerDiff + b2[Subj]*KKPerSum + b3*sex,
    b0[Subj] ~ normal(b0_mu, b0_sigma),
    b1[Subj] ~ normal(b1_mu, b1_sigma),
    b2[Subj] ~ normal(b2_mu, b2_sigma),
    b0_mu ~ normal(0,1),
    b0_sigma ~ dcauchy(0,1),
    b1_mu ~ normal(0,1),
    b1_sigma ~ dcauchy(0,1),
    b2_mu ~ normal(0,1),
    b2_sigma ~ dcauchy(0,1),
    b3 ~ normal(0,1)
  ),
  data=data,
  constraints = list(b0_sigma = "lower=0",b1_sigma = "lower=0", b2_sigma = "lower=0"),
  warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

compare(m1a,m10)
  
  
  
#########################################################################
# RT - transformed to signed response speed
#########################################################################

# Remove RT outliers
# There are some very long Rts >4000ms and some very short RTs<1ms. 
# The very long and short RTs are a problem
 dataset=data
 data = dataset[which(dataset$RT<=2000),]
 Feinting_numout=nrow(dataset)-nrow(data)
 Perout_Feinting=(Feinting_numout/nrow(dataset))*100
 print(Feinting_numout)
 print(Perout_Feinting)
 data = data[which(data$RT>=100),]
 tot_numout=nrow(dataset)-nrow(data)
 Perout_ant=((tot_numout-Feinting_numout)/nrow(dataset))*100
 print(tot_numout)
 print(Perout_ant)

rm(dataset,Feinting_numout,Perout_Feinting,tot_numout,Perout_ant)
 

# m11 regression of SRS on Logit P
# mean SRS and Logit P over all participants by condition

SummaryD=Summary(data) # function in Study1_DataPrep.R

m11 <- ulam(
  alist(
    logitP ~ normal( mu, sigma ),
    mu <- b0 + b1*SRS,
    c(b0,b1) ~ dnorm(0, 10),
    sigma ~ dunif(0 , 10)
  ) ,
  data=SummaryD, chains=chains, cores=nCores,  iter=5000, warmup=2000, sample=TRUE, log_lik=TRUE)

precis(m11,depth=2,prob=0.95)
post.m11 <- extract.samples(m11)

##################################################################
# Plot relations between SRS and Logit P
##################################################################

SRS_seq = seq(from=-1, to=1, by=0.1)
pred.data <- list(SRS=SRS_seq)
mu <- link(m11, data=pred.data)
mu.mean <- apply(mu,2,mean)
mu.HPDI <- apply(mu,2,HPDI, prob=0.95)

windows() #use X11 for mac
par(mar = c(5, 5, 2, 2))
plot(0, 0, type="n", xlab=expression(paste("Signed response speed (",italic(SRS),")")), ylab=expression(paste("Logit ", italic(P))),xlim=c(-1, 1), ylim=c(-1, 1),axes=T,bty="n", pch = 16)
lines(SRS_seq, mu.mean)
polygon(c(SRS_seq,rev(SRS_seq)),c(mu.HPDI[1,],rev(mu.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
with(SummaryD, points(SRS, logitP))

rm(mu.mean, SRS_seq, mu, mu.HPDI, pred.data)

#####################################################
#m12 SRS on Keeper, kicker, random intercepts and slopes

m12 <- ulam(
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

precis(m12,depth=3,prob=0.95)
post.m12 <- extract.samples(m12)

###########################################
# counterfactual plots SRS on keeper, kicker
##########################################
windows(height=15,width=30) #use X11 for mac
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))

# Panel A conterfactual of SRS on keeper, kicker position held constant at 0
KeeperPer_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KeeperPer) post.m12$b0_mu + post.m12$b1_mu*KeeperPer
p <- sapply(KeeperPer_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)

plot(0,0,type="n",xlab="Goalkeeper position",ylab=expression(paste("Signed response speed (",italic(SRS),")")), yaxt="n", ylim =c(-1,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6","-5","-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(-1, 1 , 0.2))
lines(KeeperPer_seq,p.mean)
polygon(c(KeeperPer_seq,rev(KeeperPer_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.0, lty=2)
segments(x0=0,y0=-1.2,x1=0,y1=1, lty=2)

SRS=with(data, tapply(SRS, KeeperPer, mean))
points(c(-4.9, -3.2, -1.6, 0, 1.6, 3.2, 4.9),SRS)

mtext(text="(% left minus right difference in goalmouth area)",side=1,line=3.5)
mtext("Kicker aligned central to the veridical goalmouth", side=3, line=1, font=2, cex = 1.5)

# Panel B conterfactual SRS on kicker, keeper position held constant at 0
Kicker_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KickerPer) post.m12$b0_mu - post.m12$b2_mu*KickerPer
p <- sapply(Kicker_seq, p.link)
p.mean <- apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)

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
mtext("Goalkeeper aligned central to the veridical goalmouth", side=3, font=2, line=1, cex=1.5)


#####################################################

#m13 SRS  on (Keeper-kicker) and (Keeper+kicker), random intercepts and slopes
m13 <- ulam(
  alist(
    SRS ~ normal( mu, sigma ),
    mu <- b0[Subj] + b1[Subj]*KKPerDiff + b2[Subj]*KKPerSum,
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

precis(m13,depth=3,prob=0.95,corr=TRUE)
post.m13 <- extract.samples(m13)


###########################################
#counterfactual plots SRS on (keeper-kicker) and (Keeper+kicker)
##########################################
windows(height=15,width=30) #use X11 for mac
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))

# Panel A conterfactual SRS on (keeperPer-kicker), (Keeper+kicker) position held constant at 0
KK_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KKPerDiff) post.m13$b0_mu + post.m13$b1_mu*KKPerDiff
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)


plot(0,0,type="n",xlab="Goalkeeper position - Kicker position",ylab=expression(paste("Signed response speed (",italic(SRS),")")), yaxt="n", ylim =c(-1,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6", "-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(-1, 1 , 0.2))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.0, lty=2)
segments(x0=0,y0=-1.2,x1=0,y1=1, lty=2)

SRS=with(data, tapply(SRS, KKPerDiff, mean))
points(c(-4.8, -1.6, 1.6, 4.8),SRS)

mtext("Goalkeeper position + Kicker position = 0", side=3, line=1, font=2, cex=1.5)

# Panel B conterfactual SRS on (Keeper+kicker), (keeper-kicker) position held constant at 0
KK_seq <- seq(from=-6, to=6, by=0.5)
p.link <- function(KKPerSum) post.m13$b0_mu + post.m13$b2_mu*KKPerSum
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI, prob=0.95)


plot(0,0,type="n",xlab="Goalkeeper position + Kicker position",ylab=expression(paste("Signed response speed (",italic(SRS),")")), yaxt="n", ylim =c(-1,1), xaxt="n", xlim=c(-6,6),bty="n")
axis(1, at=-6:6, labels=c("-6", "-5", "-4","-3","-2","-1","0","+1","+2","+3","+4","+5","+6"))
axis(2, at=seq(-1, 1 , 0.2))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.0, lty=2)
segments(x0=0,y0=-1.2,x1=0,y1=1, lty=2)

SRS=with(data, tapply(SRS, KKPerSum, mean))
points(c(-4.8, -1.6, 1.6, 4.8),SRS)

mtext("Goalkeeper position - Kicker position = 0", side=3, font=2, line=1, cex=1.5)


#######################################################################
# That's it
# Have fun modifying the script to try out different possibilities.
#######################################################################



