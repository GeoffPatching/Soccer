

SummaryBySubj <- function(d){ 
  
  #calculate Logit P for each image by subj
  for(i in min(as.numeric(d$Subj)):max(as.numeric(d$Subj))){
    temp=d[which(d$Subj==i),]
    
    
    if (exists("Im.dep")){
      table=table(temp$image,temp$Left)
      if (ncol(table)==1){logit.P=rep("-Inf",nlevels(as.factor(temp$image)))}
      
      if (ncol(table)==2){
        logit.P=log((table[,2]/rowSums(table))/(1-(table[,2]/rowSums(table))))
        PLeft=table[,2]/rowSums(table)}
         
      
      SRS=with(temp, tapply(SRS, image, mean))
      tempdata=cbind(logit.P,SRS,PLeft)
      Im.dep=rbind(Im.dep, tempdata)
      rm(tempdata,logit.P,SRS,table)
    }
    
    if (!exists("Im.dep")){
      table=table(temp$image,temp$Left)
      
      logit.P=log((table[,2]/rowSums(table))/(1-(table[,2]/rowSums(table))))
      PLeft=table[,2]/rowSums(table)
      SRS=with(temp, tapply(SRS, image, mean))
      Im.dep=cbind(logit.P,SRS,PLeft)
      rm(logit.P,SRS,table)
    }
    
    rm(temp)
  }

  Im.dep=as.data.frame(Im.dep)
  Im.dep$PLeft=as.numeric(as.character(Im.dep$PLeft))
  Im.dep$logit.P=as.numeric(as.character(Im.dep$logit.P))
  Im.dep$SRS=as.numeric(as.character(Im.dep$SRS))

  tempA <- d[!duplicated(d[,c('image')]),]
  tempA <- tempA[order(tempA$image),]
  tempB <- tempA[,c('image','KeeperPer','KickerPer')]
  n=max(as.numeric(d$Subj))
  Subj <- rep(c(1:n), each=nrow(tempB))
  tempB <- tempB[rep(seq_len(nrow(tempB)), n), ]
  tempB <- cbind(Subj, tempB)

  Im.dep<- cbind(Im.dep, tempB)
  rownames(Im.dep) <- NULL

  return(Im.dep)
}


######################################################################
Summary <- function(d){ 

  # mean SRS and Logit over all participants
  table=table(d$image,d$Left)
  tempA=d[!duplicated(d[,c('image')]),]
  tempA=tempA[order(tempA$image),]
  tempB<-tempA[,c('image','KeeperPer','KickerPer')]
  logitP=log((table[,2]/rowSums(table))/(1-(table[,2]/rowSums(table))))
  PLeft=table[,2]/rowSums(table)
  Summary.d<-cbind(tempB,logitP,PLeft)
  SRS=with(data, tapply(SRS, image, mean))
  Summary.d<-cbind(Summary.d,SRS)
  

  return(Summary.d)

}

###################################################################
# Cleans up the raw data and creates transformed variables


tidy_up=function(dataset){
  
  require(psych)
  
  #name the columns
  colnames(dataset) <- c("trials","Subj","age","finger","sex","image", "SC_keeperpos", "SC_goalpos","SC_ave", "SC_diff", "Left", "Right", "RT")
  
  #correct for MATLAB rounding error
  dataset$SC_ave[dataset$SC_ave==154.44] <- 154.43
  dataset$SC_diff[dataset$SC_diff==1.14] <- 1.13
  dataset$SC_diff[dataset$SC_diff==-1.14] <- -1.13
  
  #correct for erroroneous coding of screen co-ordinates in GP's data
  dataset$SC_keeperpos[dataset$SC_keeperpos==152.96] <- 151.60
  dataset$SC_keeperpos[dataset$SC_keeperpos==153.64] <- 152.73
  dataset$SC_keeperpos[dataset$SC_keeperpos==154.32] <- 153.87
  dataset$SC_keeperpos[dataset$SC_keeperpos==155.68] <- 156.13
  dataset$SC_keeperpos[dataset$SC_keeperpos==155.68] <- 156.13
  dataset$SC_keeperpos[dataset$SC_keeperpos==156.36] <- 157.27
  dataset$SC_keeperpos[dataset$SC_keeperpos==157.04] <- 158.40
  
  dataset$SC_goalpos[dataset$SC_goalpos==152.96] <- 151.60
  dataset$SC_goalpos[dataset$SC_goalpos==153.64] <- 152.73
  dataset$SC_goalpos[dataset$SC_goalpos==154.32] <- 153.87
  dataset$SC_goalpos[dataset$SC_goalpos==155.68] <- 156.13
  dataset$SC_goalpos[dataset$SC_goalpos==155.68] <- 156.13
  dataset$SC_goalpos[dataset$SC_goalpos==156.36] <- 157.27
  dataset$SC_goalpos[dataset$SC_goalpos==157.04] <- 158.40
  
  
  dataset$SC_ave[dataset$SC_ave==153.98] <- 153.30
  dataset$SC_ave[dataset$SC_ave==154.66] <- 154.43
  dataset$SC_ave[dataset$SC_ave==155.34] <- 155.57
  dataset$SC_ave[dataset$SC_ave==156.02] <- 156.70
  
  dataset$SC_diff[dataset$SC_diff==-0.68] <- -1.13
  dataset$SC_diff[dataset$SC_diff==-2.04] <- -3.40
  dataset$SC_diff[dataset$SC_diff==0.68] <- 1.13
  dataset$SC_diff[dataset$SC_diff==2.04] <- 3.40
  
  #check to see what we have
  dataset$Subj=as.factor(dataset$Subj)
  print(levels(as.factor(dataset$Subj)))
  print(levels(as.factor(dataset$image)))
  print(levels(as.factor(dataset$SC_goalpos)))
  print(levels(as.factor(dataset$SC_keeperpos)))
  print(levels(as.factor(dataset$SC_ave)))
  print(levels(as.factor(dataset$SC_diff)))
  print(describe(dataset$age))
  
  # #describeBy(dataset$RT,dataset$image)
  # binom.test(sum(dataset$Left),nrow(dataset))
  
  
  #recode keeper and kicker screen co-ordinates to displacements from central
  dataset$kickerpos=round(155-dataset$SC_goalpos,2)
  dataset$keeperpos=round(dataset$SC_keeperpos-155,2)
  dataset$ave=round((dataset$keeperpos+dataset$kickerpos)/2,2)
  dataset$diff=round(dataset$keeperpos-dataset$kickerpos,2)
  
  
  #correct for rounding error
  dataset$ave[dataset$ave==-0.57] <- -0.56
  dataset$ave[dataset$ave==0.57] <- 0.56
  dataset$diff[dataset$diff==1.14] <- 1.13
  dataset$diff[dataset$diff==-1.14] <- -1.13
  
  
  print(levels(as.factor(dataset$kickerpos)))
  print(levels(as.factor(dataset$keeperpos)))
  print(levels(as.factor(dataset$ave)))
  print(levels(as.factor(dataset$diff)))
  
  

  
  #difference in area 
  
  Tot_area=140*49 #physical width mm * height mm
  
  GK_Arealeft=(70-dataset$keeperpos)*49
  GK_Arearight=(70+dataset$keeperpos)*49
  
  K_Arealeft=(70-dataset$kickerpos)*49
  K_Arearight=(70+dataset$kickerpos)*49
  
  dataset$KeeperPer=(GK_Arearight-GK_Arealeft)/Tot_area
  dataset$KickerPer=(K_Arearight-K_Arealeft)/Tot_area
  
  dataset$KeeperPer=round(dataset$KeeperPer*100,1)
  dataset$KickerPer=round(dataset$KickerPer*100,1)
  
  print(levels(as.factor(dataset$KeeperPer)))
  print(levels(as.factor(dataset$KickerPer)))
  
  
  # calculate sum and difference in % difference in area
  
  dataset$KKPerSum=dataset$KeeperPer+dataset$KickerPer
  dataset$KKPerDiff=dataset$KeeperPer-dataset$KickerPer
  
  # correct for rounding error
  dataset$KKPerSum[dataset$KKPerSum==-4.9] <- -4.8
  dataset$KKPerSum[dataset$KKPerSum==4.9] <- 4.8
  dataset$KKPerDiff[dataset$KKPerDiff==4.9] <- 4.8
  dataset$KKPerDiff[dataset$KKPerDiff==-4.9] <- -4.8
  
  
  print(levels(as.factor(dataset$KKPerSum)))
  print(levels(as.factor(dataset$KKPerDiff)))
  print(levels(as.factor(dataset$KKPer_SumAve)))
  print(levels(as.factor(dataset$KKPer_DiffAve)))
  
  #recode finger and sex (male =0 , female=1)
  dataset$finger=dataset$finger-1
  dataset$sex=dataset$sex-1
  
  # Remove RT outliers
  # There are some very long Rts >4000ms and some very short RTs<1ms. 
  # The very long and short RTs are a problem
  # data = dataset[which(dataset$RT<=2000),]
  # Feinting_numout=nrow(dataset)-nrow(data)
  # Perout_Feinting=(Feinting_numout/nrow(dataset))*100
  # print(Feinting_numout)
  # print(Perout_Feinting)
  # data = data[which(data$RT>=100),]
  # tot_numout=nrow(dataset)-nrow(data)
  # Perout_ant=((tot_numout-Feinting_numout)/nrow(dataset))*100
  # print(tot_numout)
  # print(Perout_ant)
   data = dataset
  
  
  #describeBy(data$RT,data$image)
  
  #calculate SRS with RT converted to seconds - right=1 and left =-1
  data$sign=data$Left + (data$Right*-1)
  data$RTsec=data$RT/1000
  data$SRS=data$sign/data$RTsec
  # t.test(data$SRS,mu=0)
  
  return(data)
}