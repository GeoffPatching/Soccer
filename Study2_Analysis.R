# clears workspace
graphics.off()
rm(list=ls(all=TRUE)) 

# set working directory
setwd("~/Soccer")

# Import csv file as data frame
Data <- read.csv("Study2_data.csv", header=FALSE, sep=";")

# Variable Names
# Penalty - Penalty ID
# Kpos_MMI - Kickers' Position[mm] [Image]
# K_SS - Kickers' Starting Side [0= RIGHT, 1= LEFT from the goalkeepers perspective]
# GKpos_MMI - Goalkeepers' Position [Image]
# GK_dive - Goalkeeper Position Dived [0= RIGHT, 1= LEFT from the goalkeepers perspective]
# K_foot - Kickers Foot [0= RIGHT, 1= LEFT]
# KgoalSS - Kickers Goal Side Selection [0= RIGHT, 1= LEFT from the goalkeepers perspective]
# Goal - Goal ? [0=Yes, 1=No]
# Goalsize_P - Goal Size on Picture (Horizontally)
# GK_dis - Goalkeeper displacement side [0= RIGHT, 1= LEFT from the goalkeepers perspective]
colnames(Data) <- c("Penalty", "Kpos_MMI","K_SS","GKpos_MMI", "GK_dive", "K_foot", "KgoalSS", "Goal", "Goalsize_P", "GK_dis")


# Goal size official specification:
# uprights are 8 yards, or 24 feet, apart (8 yards [24 feet] = 7315.2 mm)
# Goal mouth height is 2 yards or 8 feet = 2438.4 mm   

OGS_w = 7315.2
OGS_h = 2438.4

Tot_area=OGS_w*OGS_h #physical width mm * height mm
Tot_area=Tot_area/100 # converted to cm

# Kpos - kicker position scaled relative to official goal size
# GKpos - Goalkeeper position scaled relative to official goal size
Data$Kpos = (OGS_w / Data[,"Goalsize_P"]) * Data[,"Kpos_MMI"]
Data$GKpos  = (OGS_w / Data[,"Goalsize_P"]) * Data[,"GKpos_MMI"]

# transform to kickers perspective
Data=transform(Data, Kpos_C=ifelse(K_SS==0, Kpos*1, Kpos*-1)) #1 right, -1 left, kicker perspective
Data=transform(Data, GKpos_C=ifelse(GK_dis==0, GKpos*1, GKpos*-1)) #1 right, -1 left, kicker perspective
Data$KgoalSS=1-Data$KgoalSS #0 right, 1 left, kicker perspective
Data$GK_dive=1-Data$GK_dive #0 right, 1 left, kicker perspective
Data$GK_dis=1-Data$GK_dis #0 right, 1 left, kicker perspective

# convert millimeters to centimeters
Data$Kpos_C = Data[,"Kpos_C"] / 10
Data$GKpos_C = Data[,"GKpos_C"] /10

OGS_h=OGS_h/10
OGS_w=OGS_w/10


# area to left and right of keeper and kicker
GK_Arealeft=((OGS_w/2) - Data$GKpos_C)*OGS_h
GK_Arearight=((OGS_w/2) + Data$GKpos_C)*OGS_h

K_Arealeft=((OGS_w/2) - Data$Kpos_C)*OGS_h
K_Arearight=((OGS_w/2) + Data$Kpos_C)*OGS_h

# percent change in left minus right area
Data$KeeperPer=(GK_Arealeft-GK_Arearight)/Tot_area
Data$KickerPer=(K_Arealeft-K_Arearight)/Tot_area

Data$KeeperPer=round(Data$KeeperPer*100,1)
Data$KickerPer=round(Data$KickerPer*100,1)


# sum and difference of the keepers and kickers positions
Data$KKPerSum=Data$KeeperPer+Data$KickerPer
Data$KKPerDiff=Data$KeeperPer-Data$KickerPer
Data$KKPer_SumAve=Data$KKPerSum/2
Data$KKPer_DiffAve=Data$KKPerDiff/2




#################################################
# Descriptives
##################################################
require(psych)

#some descriptive
describe(Data$KeeperPer)
describe(Data$KickerPer)

print("percentage goals")
100-sum(Data$Goal)

print("percentage Goalkeeper right dives from kicker's perspective")
sum(Data$GK_dive)

print("percentage kicker left start position from viewing perspective of the goal keeper")
sum(Data$K_SS)

print("percentage keeper left goal position from viewing perspective of the goal keeper")
sum(Data$GK_dis)

print("percent kicker position to right of the ball")
sum(  Data$KickerPer > 0  )

print("percent kicker ball kicked ball with left foot")
sum( Data$K_foot)

print("percent kicker ball kicked ball with foot opposing the start position")
sum(  Data$K_SS != Data$K_foot)
Data[which(Data$K_SS != Data$K_foot), ]

print("percent kicker left goal side selection from the kickers perspective")
sum(  Data$KgoalSS)

print("percent dive to wrong side")
sum(  Data$GK_dive != Data$KgoalSS)

print("percent kicked to goal side with greatest area")
sum(Data$GK_dis != Data$KgoalSS)

windows(height=20,width=30)
par( mar=0.5+c(5,4,0,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0))
par(mfrow=c(1,2))

brk <- seq(-12,12, by=2)
hist(Data$KeeperPer, breaks=brk, ylim=c(0,50), main="Goalkeeper position", xlab="% difference in left minus right goalmouth area from the goalkeeper's veridical midline", ylab="Frequency (out of 100 penalty shots)", xaxt="n", prob=F)
axis(1, at=c(-12,-10,-8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12), labels=c("-12","-10","-8","-6","-4","-2", "0","+2","+4","+6","+8","+10","+12"))
mtext(text="(from the kicker's perspective)",side=1,line=3.5)

brk <- seq(-140,140, by=20)
hist(Data$KickerPer, breaks=brk, ylim=c(0,50), main="Kicker position", xlab="% difference in left minus right goalmouth area from the kicker's veridical midline", ylab="Frequency (out of 100 penalty shots)", xaxt="n", prob=F)
axis(1, at=c(-150,-100,-50, 0, 50, 100, 150), labels=c("-150","-100","-50", "0","+50","+100","+150"))
mtext(text="(from the kicker's perspective)",side=1,line=3.5)


##################################################################

require(rethinking)

# leave 1 core available for other processes
nCores = parallel::detectCores() 
chains <- nCores <- nCores-1

Data$KgoalSS = as.numeric(Data$KgoalSS)

#m1 logit P on Keeper and kicker position
m1 <- ulam(
  alist(
    KgoalSS ~ bernoulli( p ),
    logit(p) <- b0 + b1*KeeperPer + b2*KickerPer,
    b0 ~ normal(0,1),
    b1 ~ normal(0,1),
    b2 ~ normal(0,1)
  ),
  data=Data, warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m1,depth=2,prob=0.95)
post.m1 <- extract.samples(m1)


#m2 logit P on (Keeper-kicker)/2 and (Keeper+kicker)/2
m2 <- ulam(
  alist(
    KgoalSS ~ bernoulli( p ),
    logit(p) <- b0 + b1*KKPer_DiffAve + b2*KKPer_SumAve,
    b0 ~ normal(0,1),
    b1 ~ normal(0,1),
    b2 ~ normal(0,1)
  ),
  data=Data, warmup=2000 , iter=5000 , chains=chains , cores=nCores, sample=TRUE, log_lik=TRUE)

precis(m2,depth=2,prob=0.95)
post.m2 <- extract.samples(m2)



###########################################
#counterfactual plots logit P on keeperPer, kicker
##########################################
windows(height=20,width=30)
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))

# Panel A conterfactual of logit P on keeperPer, kicker position held constant at 0
KeeperPer_seq <- seq(from=-12, to=12, by=0.5)
p.link <- function(KeeperPer) post.m1$b0 + post.m1$b1*KeeperPer
p <- sapply(KeeperPer_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",main="Kicker aligned central to the veridical goalmouth", xlab="Goalkeeper position",ylab="Probability of left goal shot", yaxt="n",
     ylim =c(0,1), xaxt="n", xlim=c(-12,12),bty="n")
axis(1, at=c(-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12), labels=c("-12","-10","-8", "-6","-4","-2","0","+2","+4","+6","+8","+10","+12"))
axis(2, at=seq(0, 1 , 0.1))
lines(KeeperPer_seq,p.mean)
polygon(c(KeeperPer_seq,rev(KeeperPer_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)


# Panel c conterfactual logit P on kickerPer, keeper position held constant at 0
Kicker_seq <- seq(from=-12, to=12, by=0.5)
p.link <- function(KickerPer) post.m1$b0 - post.m1$b2*KickerPer
p <- sapply(Kicker_seq, p.link)
p.mean <- apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",main="Keeper aligned central to the veridical goalmouth", xlab="Kicker position",ylab="Probability of left goal shot", yaxt="n",
    ylim =c(0,1), xaxt="n", xlim=c(-12,12),bty="n")
axis(1, at=c(-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12), labels=c("-12","-10","-8", "-6","-4","-2","0","+2","+4","+6","+8","+10","+12"))
axis(2, at=seq(0, 1 , 0.1))
lines(Kicker_seq,p.mean)
polygon(c(Kicker_seq,rev(Kicker_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)



###########################################
#counterfactual plots logit P and SRS on (keeperPer-kicker)/2 and (Keeper+kicker)/2
##########################################

windows(height=20,width=30)
par( mar=0.5+c(5,4,2,1) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
     cex.lab=1.5 )
par(mfrow=c(1,2))


# Panel A counterfactual logit P on (keeperPer-kicker)/2, (Keeper+kicker)/2 position held constant at 0
KK_seq <- seq(from=-8, to=8, by=0.5)
p.link <- function(KKPer_DiffAve) post.m2$b0 + post.m2$b1*KKPer_DiffAve
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",main="Goalkeeper position + Kicker position = 0", xlab="(Goalkeeper position - Kicker position) / 2",
     ylab="Probability of left goal side shot", yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-8,8),bty="n")
axis(1, at=c(-8,-6,-4,-2,0,2,4,6,8), labels=c("-8","-6","-4","-2","0","+2","+4","+6","+8"))
axis(2, at=seq(0, 1 , 0.1))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)



# Panel B conterfactual logit P on (Keeper+kicker)/2, (keeperPer-kicker)/2 position held constant at 0
KK_seq <- seq(from=-8, to=8, by=0.5)
p.link <- function(KKPer_SumAve) post.m2$b0 + post.m2$b2*KKPer_SumAve
p <- sapply(KK_seq, p.link)
p.mean=apply(p, 2, mean)
p.HPDI <- apply(p, 2, HPDI)
p.mean=logistic(p.mean)
p.HPDI=logistic(p.HPDI)

plot(0,0,type="n",main="Goalkeeper position - Kicker position = 0",
     xlab="(Goalkeeper position + Kicker position) / 2",ylab="Probability of left goal side shot", yaxt="n", ylim =c(0,1), xaxt="n", xlim=c(-8,8),bty="n")
axis(1, at=c(-8,-6,-4,-2,0,2,4,6,8), labels=c("-8","-6","-4","-2","0","+2","+4","+6","+8"))
axis(2, at=seq(0, 1 , 0.1))
lines(KK_seq,p.mean)
polygon(c(KK_seq,rev(KK_seq)),c(p.HPDI[1,],rev(p.HPDI[2,])), col=col.alpha("black", 0.15), border = NA)
abline(h=0.5, lty=2)
segments(x0=0,y0=-0.2,x1=0,y1=1, lty=2)





