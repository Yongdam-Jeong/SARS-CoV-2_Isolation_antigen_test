library(plyr)
library(deSolve)
library(xlsx)
library(readxl)
library(writexl)

setwd("~/Desktop/")

################################################################################################################################################
################################################################################################################################################
######################################## Generating SARS-CoV-2 infection dynamics for COVID-19 patients ########################################
################################################################################################################################################
################################################################################################################################################


#################### Setting ####################
Tmin <- 0.0
Tmax <- 140

step_size <- 1
stime <- seq(Tmin,Tmax,step_size)

sample <- 1000 ## # of samples

First <- 8 ## First test timing from (Incubation period & From symptom onset to isolation)

###### Detection limit
DL <- log10(10^4)
# DL <- log10(10^6)

###### Infectiousness threshold
IF <- log10(10^4.5)
# IF <- log10(10^5)
# IF <- log10(10^5.5)


#################### Mathematical model ####################
Covfun<-function(pars){
  beta <- as.numeric(pars[1])
  r <- as.numeric(pars[2])
  delta <- as.numeric(pars[3])
  derivs<-function(time,y,pars){
    with(as.list(c(pars,y)),{
      dTa<--beta*Ta*V
      dV<-r*Ta*V-delta*V
      
      return(list(c(dTa,dV)))
    })
  }
  y<-c(Ta=1,V=0.01)
  
  times<-c(seq(0,Tmax,step_size))
  out<-lsoda(y=y,parms=pars,times=times,func=derivs,rtol=0.00004,atol=0)
  out2<-cbind(time=out[,1],aV=((log10(out[,3]))))
  as.data.frame(out2)
}

Iteration <- 1

Prob1 <- matrix(NA,nrow=Iteration,ncol=5)
Prob2 <- matrix(NA,nrow=Iteration,ncol=5)
Prob3 <- matrix(NA,nrow=Iteration,ncol=5)
Prob4 <- matrix(NA,nrow=Iteration,ncol=5)
Prob5 <- matrix(NA,nrow=Iteration,ncol=5)

Length1 <- data.frame()
Length2 <- data.frame()
Length3 <- data.frame()
Length4 <- data.frame()
Length5 <- data.frame()

Burden1 <- data.frame()
Burden2 <- data.frame()
Burden3 <- data.frame()
Burden4 <- data.frame()
Burden5 <- data.frame()


for (l in 1:Iteration){
  
  
  #################### Parameters ####################
  pop <- read.csv("populationParameters_sym.txt", row.names = 1)
  
  popr <- rlnorm(sample, log(pop["r_pop", "value"]), pop["omega_r", "value"])  ### Log-normal dist
  popbeta <- rlnorm(sample, log(pop["beta_pop", "value"]), pop["omega_beta", "value"])
  popdelta <- rlnorm(sample, log(pop["delta_pop", "value"]), pop["omega_delta", "value"])
  fitt<-data.frame(r=popr,delta=popdelta,beta=popbeta)
  
  
  #################### Sampling 1000 COVID-19 patients ####################
  mtime<-seq(0,Tmax,step_size)
  df_total = data.frame()
  
  for( i in 1:sample ) {
    pars <- c(beta=fitt$beta[i],r=fitt$r[i],delta=fitt$delta[i])
    out <- Covfun(pars)
    
    ppp<-matrix(NA,nrow=length(mtime),ncol=1)
    for(jj in 1:length(mtime)){
      ii<-jj
      a<-mtime[ii]
      dd<-out[out$time==a,]
      ppp[ii]<-dd$aV
    }
    pp<-data.frame(ID=i, time=mtime, V=ppp)
    df_total<-rbind(df_total,pp)
    
  }
  
  
  #################### Excluding non-infectious COVID-19 patients ####################
  ex <- c()
  
  for( j in 1:sample ) {
    
    df_total2<-subset(df_total,df_total$ID==j)
    
    if (max(df_total2$V)<IF){
      ex[j]<-1
    }else {
      ex[j]<-0 
    }
    
  }
  
  ex <- data.frame(ex)
  excl <- which(ex$ex==1)
  as.numeric(excl)
  excl <- data.frame(excl)
  excl <- t(excl)
  
  df_total3 <- df_total[!(df_total$ID == excl[1]), ]
  
  for( i in 1:length(excl) ) {
    
    ii <- excl[i]
    df_total3 <- df_total3[!(df_total3$ID == ii), ]
    
  }
  
  N <- 1000 - length(excl)
  NewID <- rep(1:N, each = Tmax+1)
  df_total3$ID <- NewID
  
  
  #################### Including error for viral dynamics of COVID-19 patients ####################
  T = (Tmax+1)*N   
  TM = T+Tmax+1
  error<-rnorm(T, mean=0, sd=pop$value[9])
  
  df_totalM<-df_total3  
  df_totalM2<-df_totalM$V+error
  df_totalM$obs<-df_totalM2 ## Data frame of total patients from viral dynamics model
  
  
  ############################################################################################################
  ############################################################################################################
  ######################################## Analysis for risk & burden  #######################################
  ############################################################################################################
  ############################################################################################################
  
  nomeaning<-data.frame(ID=rep(N+1,times=Tmax+1),time=mtime,V=rep("7",times=Tmax+1),obs=rep("7",times=Tmax+1))
  df_totalM<-rbind(df_totalM,nomeaning)
  
  
  ################### Time point that VL first drops below the infectiousness threshold ###################
  cri<-c()
  
  for( i in 1:TM ) {
    if (df_totalM$V[i]<IF & df_totalM$time[i]>=First){
      cri[i]<-df_totalM$time[i]
    }else {
      cri[i]<-Tmax+1 ## meaningless
    }
  }
  
  df_totalM<-cbind(df_totalM,cri)
  df_totalM$V <- as.numeric(df_totalM$V)
  df_totalM$obs <- as.numeric(df_totalM$obs)
  
  ##############################################################################################################
  ######################################## 1 Consecutive negative result #######################################
  ##############################################################################################################
  
  pcri1<-c(); pcri2<-c(); pcri3<-c(); pcri4<-c(); pcri5<-c();
  pcri11<-c(); pcri22<-c(); pcri33<-c(); pcri44<-c(); pcri55<-c();
  
  bcri1<-c(); bcri2<-c(); bcri3<-c(); bcri4<-c(); bcri5<-c();
  bcri11<-c(); bcri22<-c(); bcri33<-c(); bcri44<-c(); bcri55<-c();
  
  lcri1<-c(); lcri2<-c(); lcri3<-c(); lcri4<-c(); lcri5<-c();
  
  #################### Interval of tests = 1 day ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First &
        df_totalM$obs[i]<DL &
        df_totalM$V[i]>IF
    ){
      pcri1[i]<-1
    }else {
      pcri1[i]<-0
    }
  }
  
  pcri1 = data.frame(pcri1)
  df_totalM<-cbind(df_totalM,pcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri1 >= 1) >=1){
      pcri11[j]<-1
    }else {
      pcri11[j]<-0
    }
  }
  
  pcri11 = data.frame(pcri11)
  prob1<-sum(pcri11==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        #& df_totalM$obs[i+1]<DL
        #& df_totalM$obs[i+2]<DL 
        #& df_totalM$obs[i+3]<DL 
        #& df_totalM$obs[i+4]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+1] 
        #& df_totalM$time[i+1]<df_totalM$time[i+2] 
        #& df_totalM$time[i+2]<df_totalM$time[i+3]
        #& df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      bcri1[i]<-df_totalM$time[i]
    }else {
      bcri1[i]<-30000 ## meaningless
    }
  }
  
  bcri1 = data.frame(bcri1)
  df_totalM<-cbind(df_totalM,bcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri11[j]<-min(df_totalM3$bcri1) - min(df_totalM3$cri)
    
  }
  
  bcri11 = data.frame(bcri11)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri11$pcri11[j]==1){
      lcri1[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri1)
    }else {
      lcri1[j]<-NaN
    }
  }
  
  lcri1 = data.frame(lcri1)
  lcri1$lcri1[lcri1$lcri1 < 0] <- 0
  lcri1 <- na.omit(lcri1)
  if (length(lcri1$lcri1) > 0) {
    lcri11 <- mean(lcri1$lcri1)
  }else {
    lcri11 <- 0
  }
  
  
  #################### Interval of tests = 2 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First & 
        df_totalM$obs[i]<DL & 
        df_totalM$V[i]>IF &
        i%%2==1
    ){
      pcri2[i]<-1
    }else {
      pcri2[i]<-0
    }
  }
  
  pcri2 = data.frame(pcri2)
  df_totalM<-cbind(df_totalM,pcri2)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri2 >= 1) >=1){
      pcri22[j]<-1
    }else {
      pcri22[j]<-0
    }
  }
  
  pcri22 = data.frame(pcri22)
  prob2<-sum(pcri22==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & i%%2==1
        #& df_totalM$obs[i+2]<DL
        #& df_totalM$obs[i+4]<DL 
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+8]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+2] 
        #& df_totalM$time[i+2]<df_totalM$time[i+4] 
        #& df_totalM$time[i+4]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      bcri2[i]<-df_totalM$time[i]
    }else {
      bcri2[i]<-30000
    }
  }
  
  bcri2 = data.frame(bcri2)
  df_totalM<-cbind(df_totalM,bcri2)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri22[j]<-min(df_totalM3$bcri2) - min(df_totalM3$cri)
    
  }
  
  bcri22 = data.frame(bcri22)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri22$pcri22[j]==1){
      lcri2[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri2)
    }else {
      lcri2[j]<-NaN
    }
  }
  
  lcri2 = data.frame(lcri2)
  lcri2$lcri2[lcri2$lcri2 < 0] <- 0
  lcri2 <- na.omit(lcri2)
  if (length(lcri2$lcri2) > 0) {
    lcri22 <- mean(lcri2$lcri2)
  }else {
    lcri22 <- 0
  }
  
  
  #################### Interval of tests = 3 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First &
        df_totalM$obs[i]<DL & 
        df_totalM$V[i]>IF &
        i%%3==1
    ){
      pcri3[i]<-1
    }else {
      pcri3[i]<-0
    }
  }
  
  pcri3 = data.frame(pcri3)
  df_totalM<-cbind(df_totalM,pcri3)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri3 >= 1) >=1){
      pcri33[j]<-1
    }else {
      pcri33[j]<-0
    }
  }
  
  pcri33 = data.frame(pcri33)
  prob3<-sum(pcri33==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%3==1
        #& df_totalM$obs[i+3]<DL
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+9]<DL 
        #& df_totalM$obs[i+12]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+3] 
        #& df_totalM$time[i+3]<df_totalM$time[i+6] 
        #& df_totalM$time[i+6]<df_totalM$time[i+9]
        #& df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      bcri3[i]<-df_totalM$time[i]
    }else {
      bcri3[i]<-30000
    }
  }
  
  bcri3 = data.frame(bcri3)
  df_totalM<-cbind(df_totalM,bcri3)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri33[j]<-min(df_totalM3$bcri3) - min(df_totalM3$cri)
    
  }
  
  bcri33 = data.frame(bcri33)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri33$pcri33[j]==1){
      lcri3[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri3)
    }else {
      lcri3[j]<-NaN
    }
  }
  
  lcri3 = data.frame(lcri3)
  lcri3$lcri3[lcri3$lcri3 < 0] <- 0
  lcri3 <- na.omit(lcri3)
  if (length(lcri3$lcri3) > 0) {
    lcri33 <- mean(lcri3$lcri3)
  }else {
    lcri33 <- 0
  }
  
  
  #################### Interval of tests = 4 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First &
        df_totalM$obs[i]<DL & 
        df_totalM$V[i]>IF &
        i%%4==1
    ){
      pcri4[i]<-1
    }else {
      pcri4[i]<-0
    }
  }
  
  pcri4 = data.frame(pcri4)
  df_totalM<-cbind(df_totalM,pcri4)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri4 >= 1) >=1){
      pcri44[j]<-1
    }else {
      pcri44[j]<-0
    }
  }
  
  pcri44 = data.frame(pcri44)
  prob4<-sum(pcri44==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & i%%4==1
        #& df_totalM$obs[i+4]<DL
        #& df_totalM$obs[i+8]<DL 
        #& df_totalM$obs[i+12]<DL 
        #& df_totalM$obs[i+16]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+4] 
        #& df_totalM$time[i+4]<df_totalM$time[i+8] 
        #& df_totalM$time[i+8]<df_totalM$time[i+12]
        #& df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      bcri4[i]<-df_totalM$time[i]
    }else {
      bcri4[i]<-30000
    }
  }
  
  bcri4 = data.frame(bcri4)
  df_totalM<-cbind(df_totalM,bcri4)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri44[j]<-min(df_totalM3$bcri4) - min(df_totalM3$cri)
    
  }
  
  bcri44 = data.frame(bcri44)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri44$pcri44[j]==1){
      lcri4[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri4)
    }else {
      lcri4[j]<-NaN
    }
  }
  
  lcri4 = data.frame(lcri4)
  lcri4$lcri4[lcri4$lcri4 < 0] <- 0
  lcri4 <- na.omit(lcri4)
  if (length(lcri4$lcri4) > 0) {
    lcri44 <- mean(lcri4$lcri4)
  }else {
    lcri44 <- 0
  }
  
  
  #################### Interval of tests = 5 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First &
        df_totalM$obs[i]<DL & 
        df_totalM$V[i]>IF &
        i%%5==1
    ){
      pcri5[i]<-1
    }else {
      pcri5[i]<-0
    }
  }
  
  pcri5 = data.frame(pcri5)
  df_totalM<-cbind(df_totalM,pcri5)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri5 >= 1) >=1){
      pcri55[j]<-1
    }else {
      pcri55[j]<-0
    }
  }
  
  pcri55 = data.frame(pcri55)
  prob5<-sum(pcri55==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & i%%5==1
        #& df_totalM$obs[i+5]<DL
        #& df_totalM$obs[i+10]<DL 
        #& df_totalM$obs[i+15]<DL 
        #& df_totalM$obs[i+20]<DL 
        #& df_totalM$time[i]<df_totalM$time[i+5] 
        #& df_totalM$time[i+5]<df_totalM$time[i+10] 
        #& df_totalM$time[i+10]<df_totalM$time[i+15]
        #& df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      bcri5[i]<-df_totalM$time[i]
    }else {
      bcri5[i]<-30000
    }
  }
  
  bcri5 = data.frame(bcri5)
  df_totalM<-cbind(df_totalM,bcri5)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri55[j]<-min(df_totalM3$bcri5) - min(df_totalM3$cri)
    
  }
  
  bcri55 = data.frame(bcri55)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri55$pcri55[j]==1){
      lcri5[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri5)
    }else {
      lcri5[j]<-NaN
    }
  }
  
  lcri5 = data.frame(lcri5)
  lcri5$lcri5[lcri5$lcri5 < 0] <- 0
  lcri5 <- na.omit(lcri5)
  if (length(lcri5$lcri5) > 0) {
    lcri55 <- mean(lcri5$lcri5)
  }else {
    lcri55 <- 0
  }
  
  
  #################### Result ####################
  Prob1[l,] <- c(prob1,prob2,prob3,prob4,prob5)
  
  B1 <- data.frame(bcri11,bcri22,bcri33,bcri44,bcri55)
  Burden1 <- rbind(Burden1,B1)
  
  L1 <- data.frame(lcri11,lcri22,lcri33,lcri44,lcri55)
  Length1 <- rbind(Length1,L1)
  
  
  df_totalM<-subset(df_totalM,select=-c(pcri1,pcri2,pcri3,pcri4,pcri5,bcri1,bcri2,bcri3,bcri4,bcri5))
  pcri1<-c(); pcri2<-c(); pcri3<-c(); pcri4<-c(); pcri5<-c();
  pcri11<-c(); pcri22<-c(); pcri33<-c(); pcri44<-c(); pcri55<-c();
  bcri1<-c(); bcri2<-c(); bcri3<-c(); bcri4<-c(); bcri5<-c();
  bcri11<-c(); bcri22<-c(); bcri33<-c(); bcri44<-c(); bcri55<-c();
  lcri1<-c(); lcri2<-c(); lcri3<-c(); lcri4<-c(); lcri5<-c();
  rm(df_totalM3)
  
  
  
  
  
  ###############################################################################################################
  ######################################## 2 Consecutive negative results #######################################
  ###############################################################################################################
  
  # nomeaning<-data.frame(ID=rep(N+1,times=Tmax+1),time=mtime,V=rep("7",times=Tmax+1),obs=rep("7",times=Tmax+1))
  # df_totalM<-rbind(df_totalM,nomeaning)
  
  #################### Interval of tests = 1 day ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First 
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL 
        #& df_totalM$obs[i+2]<DL 
        #& df_totalM$obs[i+3]<DL 
        #& df_totalM$obs[i+4]<DL
        #& df_totalM$V[i]>DL
        & df_totalM$V[i+1]>IF
        #& df_totalM$V[i+2]>DL 
        #& df_totalM$V[i+3]>DL
        #& df_totalM$V[i+4]>DL 
        & df_totalM$time[i]<df_totalM$time[i+1] 
        #& df_totalM$time[i+1]<df_totalM$time[i+2] 
        #& df_totalM$time[i+2]<df_totalM$time[i+3]
        #& df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      pcri1[i]<-1
    }else {
      pcri1[i]<-0
    }
  }
  
  pcri1 = data.frame(pcri1)
  df_totalM<-cbind(df_totalM,pcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri1 >= 1) >=1){
      pcri11[j]<-1
    }else {
      pcri11[j]<-0
    }
  }
  
  pcri11 = data.frame(pcri11)
  prob1<-sum(pcri11==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL
        #& df_totalM$obs[i+2]<DL 
        #& df_totalM$obs[i+3]<DL 
        #& df_totalM$obs[i+4]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1] 
        #& df_totalM$time[i+1]<df_totalM$time[i+2] 
        #& df_totalM$time[i+2]<df_totalM$time[i+3]
        #& df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      bcri1[i]<-df_totalM$time[i+1]
    }else {
      bcri1[i]<-30000
    }
  }
  
  bcri1 = data.frame(bcri1)
  df_totalM<-cbind(df_totalM,bcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri11[j]<-min(df_totalM3$bcri1) - min(df_totalM3$cri)
    
  }
  
  bcri11 = data.frame(bcri11)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri11$pcri11[j]==1){
      lcri1[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri1)
    }else {
      lcri1[j]<-NaN
    }
  }
  
  lcri1 = data.frame(lcri1)
  lcri1$lcri1[lcri1$lcri1 < 0] <- 0
  lcri1 <- na.omit(lcri1)
  if (length(lcri1$lcri1) > 0) {
    lcri11 <- mean(lcri1$lcri1)
  }else {
    lcri11 <- 0
  }
  
  
  #################### Interval of tests = 2 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+2]<DL 
        #& df_totalM$obs[i+4]<DL 
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+8]<DL 
        & df_totalM$V[i+2]>IF
        #& df_totalM$V[i+4]>DL
        #& df_totalM$V[i+6]>DL
        #& df_totalM$V[i+8]>DL
        & df_totalM$time[i]<df_totalM$time[i+2] 
        #& df_totalM$time[i+2]<df_totalM$time[i+4] 
        #& df_totalM$time[i+4]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      pcri2[i]<-1
    }else {
      pcri2[i]<-0
    }
  }
  
  pcri2 = data.frame(pcri2)
  df_totalM<-cbind(df_totalM,pcri2)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri2 >= 1) >=1){
      pcri22[j]<-1
    }else {
      pcri22[j]<-0
    }
  }
  
  pcri22 = data.frame(pcri22)
  prob2<-sum(pcri22==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+2]<DL
        #& df_totalM$obs[i+4]<DL 
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+8]<DL 
        & df_totalM$time[i]<df_totalM$time[i+2] 
        #& df_totalM$time[i+2]<df_totalM$time[i+4] 
        #& df_totalM$time[i+4]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      bcri2[i]<-df_totalM$time[i+2]
    }else {
      bcri2[i]<-30000
    }
  }
  
  bcri2 = data.frame(bcri2)
  df_totalM<-cbind(df_totalM,bcri2)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri22[j]<- min(df_totalM3$bcri2) - min(df_totalM3$cri)
    
  }
  
  bcri22 = data.frame(bcri22)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri22$pcri22[j]==1){
      lcri2[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri2)
    }else {
      lcri2[j]<-NaN
    }
  }
  
  lcri2 = data.frame(lcri2)
  lcri2$lcri2[lcri2$lcri2 < 0] <- 0
  lcri2 <- na.omit(lcri2)
  if (length(lcri2$lcri2) > 0) {
    lcri22 <- mean(lcri2$lcri2)
  }else {
    lcri22 <- 0
  }
  
  
  #################### Interval of tests = 3 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+3]<DL 
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+9]<DL 
        #& df_totalM$obs[i+12]<DL
        & df_totalM$V[i+3]>IF
        #& df_totalM$V[i+6]>DL
        #& df_totalM$V[i+9]>DL
        #& df_totalM$V[i+12]>DL
        & df_totalM$time[i]<df_totalM$time[i+3] 
        #& df_totalM$time[i+3]<df_totalM$time[i+6] 
        #& df_totalM$time[i+6]<df_totalM$time[i+9]
        #& df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      pcri3[i]<-1
    }else {
      pcri3[i]<-0
    }
  }
  
  pcri3 = data.frame(pcri3)
  df_totalM<-cbind(df_totalM,pcri3)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri3 >= 1) >=1){
      pcri33[j]<-1
    }else {
      pcri33[j]<-0
    }
  }
  
  pcri33 = data.frame(pcri33)
  prob3<-sum(pcri33==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+3]<DL
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+9]<DL 
        #& df_totalM$obs[i+12]<DL 
        & df_totalM$time[i]<df_totalM$time[i+3] 
        #& df_totalM$time[i+3]<df_totalM$time[i+6] 
        #& df_totalM$time[i+6]<df_totalM$time[i+9]
        #& df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      bcri3[i]<-df_totalM$time[i+3]
    }else {
      bcri3[i]<-30000
    }
  }
  
  bcri3 = data.frame(bcri3)
  df_totalM<-cbind(df_totalM,bcri3)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri33[j]<-min(df_totalM3$bcri3) - min(df_totalM3$cri)
    
  }
  
  bcri33 = data.frame(bcri33)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri33$pcri33[j]==1){
      lcri3[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri3)
    }else {
      lcri3[j]<-NaN
    }
  }
  
  lcri3 = data.frame(lcri3)
  lcri3$lcri3[lcri3$lcri3 < 0] <- 0
  lcri3 <- na.omit(lcri3)
  if (length(lcri3$lcri3) > 0) {
    lcri33 <- mean(lcri3$lcri3)
  }else {
    lcri33 <- 0
  }
  
  
  #################### Interval of tests = 4 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+4]<DL 
        #& df_totalM$obs[i+8]<DL 
        #& df_totalM$obs[i+12]<DL 
        #& df_totalM$obs[i+16]<DL
        & df_totalM$V[i+4]>IF
        #& df_totalM$V[i+8]>DL
        #& df_totalM$V[i+12]>DL
        #& df_totalM$V[i+16]>DL
        & df_totalM$time[i]<df_totalM$time[i+4] 
        #& df_totalM$time[i+4]<df_totalM$time[i+8] 
        #& df_totalM$time[i+8]<df_totalM$time[i+12]
        #& df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      pcri4[i]<-1
    }else {
      pcri4[i]<-0
    }
  }
  
  pcri4 = data.frame(pcri4)
  df_totalM<-cbind(df_totalM,pcri4)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri4 >= 1) >=1){
      pcri44[j]<-1
    }else {
      pcri44[j]<-0
    }
  }
  
  pcri44 = data.frame(pcri44)
  prob4<-sum(pcri44==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+4]<DL
        #& df_totalM$obs[i+8]<DL 
        #& df_totalM$obs[i+12]<DL 
        #& df_totalM$obs[i+16]<DL 
        & df_totalM$time[i]<df_totalM$time[i+4] 
        #& df_totalM$time[i+4]<df_totalM$time[i+8] 
        #& df_totalM$time[i+8]<df_totalM$time[i+12]
        #& df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      bcri4[i]<-df_totalM$time[i+4]
    }else {
      bcri4[i]<-30000
    }
  }
  
  bcri4 = data.frame(bcri4)
  df_totalM<-cbind(df_totalM,bcri4)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri44[j]<-min(df_totalM3$bcri4) - min(df_totalM3$cri)
    
  }
  
  bcri44 = data.frame(bcri44)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri44$pcri44[j]==1){
      lcri4[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri4)
    }else {
      lcri4[j]<-NaN
    }
  }
  
  lcri4 = data.frame(lcri4)
  lcri4$lcri4[lcri4$lcri4 < 0] <- 0
  lcri4 <- na.omit(lcri4)
  if (length(lcri4$lcri4) > 0) {
    lcri44 <- mean(lcri4$lcri4)
  }else {
    lcri44 <- 0
  }
  
  
  #################### Interval of tests = 5 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL 
        #& df_totalM$obs[i+10]<DL 
        #& df_totalM$obs[i+15]<DL 
        #& df_totalM$obs[i+20]<DL
        & df_totalM$V[i+5]>IF
        #& df_totalM$V[i+10]>DL
        #& df_totalM$V[i+15]>DL
        #& df_totalM$V[i+20]>DL
        & df_totalM$time[i]<df_totalM$time[i+5] 
        #& df_totalM$time[i+5]<df_totalM$time[i+10] 
        #& df_totalM$time[i+10]<df_totalM$time[i+15]
        #& df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      pcri5[i]<-1
    }else {
      pcri5[i]<-0
    }
  }
  
  pcri5 = data.frame(pcri5)
  df_totalM<-cbind(df_totalM,pcri5)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri5 >= 1) >=1){
      pcri55[j]<-1
    }else {
      pcri55[j]<-0
    }
  }
  
  pcri55 = data.frame(pcri55)
  prob5<-sum(pcri55==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL
        #& df_totalM$obs[i+10]<DL 
        #& df_totalM$obs[i+15]<DL 
        #& df_totalM$obs[i+20]<DL 
        & df_totalM$time[i]<df_totalM$time[i+5] 
        #& df_totalM$time[i+5]<df_totalM$time[i+10] 
        #& df_totalM$time[i+10]<df_totalM$time[i+15]
        #& df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      bcri5[i]<-df_totalM$time[i+5]
    }else {
      bcri5[i]<-30000
    }
  }
  
  bcri5 = data.frame(bcri5)
  df_totalM<-cbind(df_totalM,bcri5)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri55[j]<-min(df_totalM3$bcri5) - min(df_totalM3$cri)
    
  }
  
  bcri55 = data.frame(bcri55)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri55$pcri55[j]==1){
      lcri5[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri5)
    }else {
      lcri5[j]<-NaN
    }
  }
  
  lcri5 = data.frame(lcri5)
  lcri5$lcri5[lcri5$lcri5 < 0] <- 0
  lcri5 <- na.omit(lcri5)
  if (length(lcri5$lcri5) > 0) {
    lcri55 <- mean(lcri5$lcri5)
  }else {
    lcri55 <- 0
  }
  
  
  #################### Result ####################
  Prob2[l,] <- c(prob1,prob2,prob3,prob4,prob5)
  
  B2 <- data.frame(bcri11,bcri22,bcri33,bcri44,bcri55)
  Burden2 <- rbind(Burden2,B2)
  
  L2 <- data.frame(lcri11,lcri22,lcri33,lcri44,lcri55)
  Length2 <- rbind(Length2,L2)
  
  
  df_totalM<-subset(df_totalM,select=-c(pcri1,pcri2,pcri3,pcri4,pcri5,bcri1,bcri2,bcri3,bcri4,bcri5))
  pcri1<-c(); pcri2<-c(); pcri3<-c(); pcri4<-c(); pcri5<-c();
  pcri11<-c(); pcri22<-c(); pcri33<-c(); pcri44<-c(); pcri55<-c();
  bcri1<-c(); bcri2<-c(); bcri3<-c(); bcri4<-c(); bcri5<-c();
  bcri11<-c(); bcri22<-c(); bcri33<-c(); bcri44<-c(); bcri55<-c();
  lcri1<-c(); lcri2<-c(); lcri3<-c(); lcri4<-c(); lcri5<-c();
  rm(df_totalM3)
  
  
  
  
  ###############################################################################################################
  ######################################## 3 Consecutive negative results #######################################
  ###############################################################################################################
  
  #################### Interval of tests = 1 day ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL 
        & df_totalM$obs[i+2]<DL 
        #& df_totalM$obs[i+3]<DL 
        #& df_totalM$obs[i+4]<DL
        #& df_totalM$V[i]>DL
        #& df_totalM$V[i+1]>DL 
        & df_totalM$V[i+2]>IF 
        #& df_totalM$V[i+3]>DL
        #& df_totalM$V[i+4]>DL 
        & df_totalM$time[i]<df_totalM$time[i+1] 
        & df_totalM$time[i+1]<df_totalM$time[i+2] 
        #& df_totalM$time[i+2]<df_totalM$time[i+3]
        #& df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      pcri1[i]<-1
    }else {
      pcri1[i]<-0
    }
  }
  
  pcri1 = data.frame(pcri1)
  df_totalM<-cbind(df_totalM,pcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri1 >= 1) >=1){
      pcri11[j]<-1
    }else {
      pcri11[j]<-0
    }
  }
  
  pcri11 = data.frame(pcri11)
  prob1<-sum(pcri11==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL
        & df_totalM$obs[i+2]<DL
        #& df_totalM$obs[i+3]<DL 
        #& df_totalM$obs[i+4]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1] 
        & df_totalM$time[i+1]<df_totalM$time[i+2]
        #& df_totalM$time[i+2]<df_totalM$time[i+3]
        #& df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      bcri1[i]<-df_totalM$time[i+2]
    }else {
      bcri1[i]<-30000
    }
  }
  
  bcri1 = data.frame(bcri1)
  df_totalM<-cbind(df_totalM,bcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri11[j]<-min(df_totalM3$bcri1) - min(df_totalM3$cri)
    
  }
  
  bcri11 = data.frame(bcri11)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri11$pcri11[j]==1){
      lcri1[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri1)
    }else {
      lcri1[j]<-NaN
    }
  }
  
  lcri1 = data.frame(lcri1)
  lcri1$lcri1[lcri1$lcri1 < 0] <- 0
  lcri1 <- na.omit(lcri1)
  if (length(lcri1$lcri1) > 0) {
    lcri11 <- mean(lcri1$lcri1)
  }else {
    lcri11 <- 0
  }
  
  
  #################### Interval of tests = 2 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+2]<DL 
        & df_totalM$obs[i+4]<DL 
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+8]<DL 
        #& df_totalM$V[i+2]>DL
        & df_totalM$V[i+4]>IF
        #& df_totalM$V[i+6]>DL
        #& df_totalM$V[i+8]>DL
        & df_totalM$time[i]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+4] 
        #& df_totalM$time[i+4]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      pcri2[i]<-1
    }else {
      pcri2[i]<-0
    }
  }
  
  pcri2 = data.frame(pcri2)
  df_totalM<-cbind(df_totalM,pcri2)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri2 >= 1) >=1){
      pcri22[j]<-1
    }else {
      pcri22[j]<-0
    }
  }
  
  pcri22 = data.frame(pcri22)
  prob2<-sum(pcri22==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+2]<DL
        & df_totalM$obs[i+4]<DL
        #& df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+8]<DL 
        & df_totalM$time[i]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+4]
        #& df_totalM$time[i+4]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      bcri2[i]<-df_totalM$time[i+4]
    }else {
      bcri2[i]<-30000
    }
  }
  
  bcri2 = data.frame(bcri2)
  df_totalM<-cbind(df_totalM,bcri2)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri22[j]<- min(df_totalM3$bcri2) - min(df_totalM3$cri)
    
  }
  
  bcri22 = data.frame(bcri22)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri22$pcri22[j]==1){
      lcri2[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri2)
    }else {
      lcri2[j]<-NaN
    }
  }
  
  lcri2 = data.frame(lcri2)
  lcri2$lcri2[lcri2$lcri2 < 0] <- 0
  lcri2 <- na.omit(lcri2)
  if (length(lcri2$lcri2) > 0) {
    lcri22 <- mean(lcri2$lcri2)
  }else {
    lcri22 <- 0
  }
  
  
  #################### Interval of tests = 3 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+3]<DL 
        & df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+9]<DL 
        #& df_totalM$obs[i+12]<DL
        #& df_totalM$V[i+3]>DL
        & df_totalM$V[i+6]>IF
        #& df_totalM$V[i+9]>DL
        #& df_totalM$V[i+12]>DL
        & df_totalM$time[i]<df_totalM$time[i+3] 
        & df_totalM$time[i+3]<df_totalM$time[i+6] 
        #& df_totalM$time[i+6]<df_totalM$time[i+9]
        #& df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      pcri3[i]<-1
    }else {
      pcri3[i]<-0
    }
  }
  
  pcri3 = data.frame(pcri3)
  df_totalM<-cbind(df_totalM,pcri3)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri3 >= 1) >=1){
      pcri33[j]<-1
    }else {
      pcri33[j]<-0
    }
  }
  
  pcri33 = data.frame(pcri33)
  prob3<-sum(pcri33==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+3]<DL
        & df_totalM$obs[i+6]<DL
        #& df_totalM$obs[i+9]<DL 
        #& df_totalM$obs[i+12]<DL 
        & df_totalM$time[i]<df_totalM$time[i+3] 
        & df_totalM$time[i+3]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+9]
        #& df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      bcri3[i]<-df_totalM$time[i+6]
    }else {
      bcri3[i]<-30000
    }
  }
  
  bcri3 = data.frame(bcri3)
  df_totalM<-cbind(df_totalM,bcri3)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri33[j]<-min(df_totalM3$bcri3) - min(df_totalM3$cri)
    
  }
  
  bcri33 = data.frame(bcri33)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri33$pcri33[j]==1){
      lcri3[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri3)
    }else {
      lcri3[j]<-NaN
    }
  }
  
  lcri3 = data.frame(lcri3)
  lcri3$lcri3[lcri3$lcri3 < 0] <- 0
  lcri3 <- na.omit(lcri3)
  if (length(lcri3$lcri3) > 0) {
    lcri33 <- mean(lcri3$lcri3)
  }else {
    lcri33 <- 0
  }
  
  
  #################### Interval of tests = 4 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+4]<DL 
        & df_totalM$obs[i+8]<DL 
        #& df_totalM$obs[i+12]<DL 
        #& df_totalM$obs[i+16]<DL
        #& df_totalM$V[i+4]>DL
        & df_totalM$V[i+8]>IF
        #& df_totalM$V[i+12]>DL
        #& df_totalM$V[i+16]>DL
        & df_totalM$time[i]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+8] 
        #& df_totalM$time[i+8]<df_totalM$time[i+12]
        #& df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      pcri4[i]<-1
    }else {
      pcri4[i]<-0
    }
  }
  
  pcri4 = data.frame(pcri4)
  df_totalM<-cbind(df_totalM,pcri4)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri4 >= 1) >=1){
      pcri44[j]<-1
    }else {
      pcri44[j]<-0
    }
  }
  
  pcri44 = data.frame(pcri44)
  prob4<-sum(pcri44==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+4]<DL
        & df_totalM$obs[i+8]<DL
        #& df_totalM$obs[i+12]<DL 
        #& df_totalM$obs[i+16]<DL 
        & df_totalM$time[i]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+8]
        #& df_totalM$time[i+8]<df_totalM$time[i+12]
        #& df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      bcri4[i]<-df_totalM$time[i+8]
    }else {
      bcri4[i]<-30000
    }
  }
  
  bcri4 = data.frame(bcri4)
  df_totalM<-cbind(df_totalM,bcri4)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri44[j]<-min(df_totalM3$bcri4) - min(df_totalM3$cri)
    
  }
  
  bcri44 = data.frame(bcri44)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri44$pcri44[j]==1){
      lcri4[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri4)
    }else {
      lcri4[j]<-NaN
    }
  }
  
  lcri4 = data.frame(lcri4)
  lcri4$lcri4[lcri4$lcri4 < 0] <- 0
  lcri4 <- na.omit(lcri4)
  if (length(lcri4$lcri4) > 0) {
    lcri44 <- mean(lcri4$lcri4)
  }else {
    lcri44 <- 0
  }
  
  
  #################### Interval of tests = 5 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL 
        & df_totalM$obs[i+10]<DL 
        #& df_totalM$obs[i+15]<DL 
        #& df_totalM$obs[i+20]<DL
        #& df_totalM$V[i+5]>DL
        & df_totalM$V[i+10]>IF
        #& df_totalM$V[i+15]>DL
        #& df_totalM$V[i+20]>DL
        & df_totalM$time[i]<df_totalM$time[i+5] 
        & df_totalM$time[i+5]<df_totalM$time[i+10] 
        #& df_totalM$time[i+10]<df_totalM$time[i+15]
        #& df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      pcri5[i]<-1
    }else {
      pcri5[i]<-0
    }
  }
  
  pcri5 = data.frame(pcri5)
  df_totalM<-cbind(df_totalM,pcri5)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri5 >= 1) >=1){
      pcri55[j]<-1
    }else {
      pcri55[j]<-0
    }
  }
  
  pcri55 = data.frame(pcri55)
  prob5<-sum(pcri55==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL
        & df_totalM$obs[i+10]<DL
        #& df_totalM$obs[i+15]<DL 
        #& df_totalM$obs[i+20]<DL 
        & df_totalM$time[i]<df_totalM$time[i+5] 
        & df_totalM$time[i+5]<df_totalM$time[i+10]
        #& df_totalM$time[i+10]<df_totalM$time[i+15]
        #& df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      bcri5[i]<-df_totalM$time[i+10]
    }else {
      bcri5[i]<-30000
    }
  }
  
  bcri5 = data.frame(bcri5)
  df_totalM<-cbind(df_totalM,bcri5)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri55[j]<-min(df_totalM3$bcri5) - min(df_totalM3$cri)
    
  }
  
  bcri55 = data.frame(bcri55)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri55$pcri55[j]==1){
      lcri5[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri5)
    }else {
      lcri5[j]<-NaN
    }
  }
  
  lcri5 = data.frame(lcri5)
  lcri5$lcri5[lcri5$lcri5 < 0] <- 0
  lcri5 <- na.omit(lcri5)
  if (length(lcri5$lcri5) > 0) {
    lcri55 <- mean(lcri5$lcri5)
  }else {
    lcri55 <- 0
  }
  
  
  #################### Result ####################
  Prob3[l,] <- c(prob1,prob2,prob3,prob4,prob5)
  
  B3 <- data.frame(bcri11,bcri22,bcri33,bcri44,bcri55)
  Burden3 <- rbind(Burden3,B3)
  
  L3 <- data.frame(lcri11,lcri22,lcri33,lcri44,lcri55)
  Length3 <- rbind(Length3,L3)
  
  
  df_totalM<-subset(df_totalM,select=-c(pcri1,pcri2,pcri3,pcri4,pcri5,bcri1,bcri2,bcri3,bcri4,bcri5))
  pcri1<-c(); pcri2<-c(); pcri3<-c(); pcri4<-c(); pcri5<-c();
  pcri11<-c(); pcri22<-c(); pcri33<-c(); pcri44<-c(); pcri55<-c();
  bcri1<-c(); bcri2<-c(); bcri3<-c(); bcri4<-c(); bcri5<-c();
  bcri11<-c(); bcri22<-c(); bcri33<-c(); bcri44<-c(); bcri55<-c();
  lcri1<-c(); lcri2<-c(); lcri3<-c(); lcri4<-c(); lcri5<-c();
  rm(df_totalM3)
  
  
  
  
  ###############################################################################################################
  ######################################## 4 Consecutive negative results #######################################
  ###############################################################################################################
  
  #################### Interval of tests = 1 day ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL 
        & df_totalM$obs[i+2]<DL 
        & df_totalM$obs[i+3]<DL 
        #& df_totalM$obs[i+4]<DL
        #& df_totalM$V[i]>DL
        #& df_totalM$V[i+1]>DL 
        #& df_totalM$V[i+2]>DL 
        & df_totalM$V[i+3]>IF
        #& df_totalM$V[i+4]>DL 
        & df_totalM$time[i]<df_totalM$time[i+1] 
        & df_totalM$time[i+1]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+3]
        #& df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      pcri1[i]<-1
    }else {
      pcri1[i]<-0
    }
  }
  
  pcri1 = data.frame(pcri1)
  df_totalM<-cbind(df_totalM,pcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri1 >= 1) >=1){
      pcri11[j]<-1
    }else {
      pcri11[j]<-0
    }
  }
  
  pcri11 = data.frame(pcri11)
  prob1<-sum(pcri11==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL
        & df_totalM$obs[i+2]<DL
        & df_totalM$obs[i+3]<DL
        #& df_totalM$obs[i+4]<DL 
        & df_totalM$time[i]<df_totalM$time[i+1] 
        & df_totalM$time[i+1]<df_totalM$time[i+2]
        & df_totalM$time[i+2]<df_totalM$time[i+3]
        #& df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      bcri1[i]<-df_totalM$time[i+3]
    }else {
      bcri1[i]<-30000
    }
  }
  
  bcri1 = data.frame(bcri1)
  df_totalM<-cbind(df_totalM,bcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri11[j]<-min(df_totalM3$bcri1) - min(df_totalM3$cri)
    
  }
  
  bcri11 = data.frame(bcri11)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri11$pcri11[j]==1){
      lcri1[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri1)
    }else {
      lcri1[j]<-NaN
    }
  }
  
  lcri1 = data.frame(lcri1)
  lcri1$lcri1[lcri1$lcri1 < 0] <- 0
  lcri1 <- na.omit(lcri1)
  if (length(lcri1$lcri1) > 0) {
    lcri11 <- mean(lcri1$lcri1)
  }else {
    lcri11 <- 0
  }
  
  
  #################### Interval of tests = 2 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+2]<DL 
        & df_totalM$obs[i+4]<DL 
        & df_totalM$obs[i+6]<DL 
        #& df_totalM$obs[i+8]<DL 
        #& df_totalM$V[i+2]>DL
        #& df_totalM$V[i+4]>DL
        & df_totalM$V[i+6]>IF
        #& df_totalM$V[i+8]>DL
        & df_totalM$time[i]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      pcri2[i]<-1
    }else {
      pcri2[i]<-0
    }
  }
  
  pcri2 = data.frame(pcri2)
  df_totalM<-cbind(df_totalM,pcri2)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri2 >= 1) >=1){
      pcri22[j]<-1
    }else {
      pcri22[j]<-0
    }
  }
  
  pcri22 = data.frame(pcri22)
  prob2<-sum(pcri22==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+2]<DL
        & df_totalM$obs[i+4]<DL
        & df_totalM$obs[i+6]<DL
        #& df_totalM$obs[i+8]<DL 
        & df_totalM$time[i]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+4]
        & df_totalM$time[i+4]<df_totalM$time[i+6]
        #& df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      bcri2[i]<-df_totalM$time[i+6]
    }else {
      bcri2[i]<-30000
    }
  }
  
  bcri2 = data.frame(bcri2)
  df_totalM<-cbind(df_totalM,bcri2)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri22[j]<- min(df_totalM3$bcri2) - min(df_totalM3$cri)
    
  }
  
  bcri22 = data.frame(bcri22)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri22$pcri22[j]==1){
      lcri2[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri2)
    }else {
      lcri2[j]<-NaN
    }
  }
  
  lcri2 = data.frame(lcri2)
  lcri2$lcri2[lcri2$lcri2 < 0] <- 0
  lcri2 <- na.omit(lcri2)
  if (length(lcri2$lcri2) > 0) {
    lcri22 <- mean(lcri2$lcri2)
  }else {
    lcri22 <- 0
  }
  
  
  #################### Interval of tests = 3 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+3]<DL 
        & df_totalM$obs[i+6]<DL 
        & df_totalM$obs[i+9]<DL 
        #& df_totalM$obs[i+12]<DL
        #& df_totalM$V[i+3]>DL
        #& df_totalM$V[i+6]>DL
        & df_totalM$V[i+9]>IF
        #& df_totalM$V[i+12]>DL
        & df_totalM$time[i]<df_totalM$time[i+3] 
        & df_totalM$time[i+3]<df_totalM$time[i+6] 
        & df_totalM$time[i+6]<df_totalM$time[i+9]
        #& df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      pcri3[i]<-1
    }else {
      pcri3[i]<-0
    }
  }
  
  pcri3 = data.frame(pcri3)
  df_totalM<-cbind(df_totalM,pcri3)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri3 >= 1) >=1){
      pcri33[j]<-1
    }else {
      pcri33[j]<-0
    }
  }
  
  pcri33 = data.frame(pcri33)
  prob3<-sum(pcri33==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+3]<DL
        & df_totalM$obs[i+6]<DL
        & df_totalM$obs[i+9]<DL
        #& df_totalM$obs[i+12]<DL 
        & df_totalM$time[i]<df_totalM$time[i+3] 
        & df_totalM$time[i+3]<df_totalM$time[i+6]
        & df_totalM$time[i+6]<df_totalM$time[i+9]
        #& df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      bcri3[i]<-df_totalM$time[i+9]
    }else {
      bcri3[i]<-30000
    }
  }
  
  bcri3 = data.frame(bcri3)
  df_totalM<-cbind(df_totalM,bcri3)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri33[j]<-min(df_totalM3$bcri3) - min(df_totalM3$cri)
    
  }
  
  bcri33 = data.frame(bcri33)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri33$pcri33[j]==1){
      lcri3[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri3)
    }else {
      lcri3[j]<-NaN
    }
  }
  
  lcri3 = data.frame(lcri3)
  lcri3$lcri3[lcri3$lcri3 < 0] <- 0
  lcri3 <- na.omit(lcri3)
  if (length(lcri3$lcri3) > 0) {
    lcri33 <- mean(lcri3$lcri3)
  }else {
    lcri33 <- 0
  }
  
  
  #################### Interval of tests = 4 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+4]<DL 
        & df_totalM$obs[i+8]<DL 
        & df_totalM$obs[i+12]<DL 
        #& df_totalM$obs[i+16]<DL
        #& df_totalM$V[i+4]>DL
        #& df_totalM$V[i+8]>DL
        & df_totalM$V[i+12]>IF
        #& df_totalM$V[i+16]>DL
        & df_totalM$time[i]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+8] 
        & df_totalM$time[i+8]<df_totalM$time[i+12]
        #& df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      pcri4[i]<-1
    }else {
      pcri4[i]<-0
    }
  }
  
  pcri4 = data.frame(pcri4)
  df_totalM<-cbind(df_totalM,pcri4)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri4 >= 1) >=1){
      pcri44[j]<-1
    }else {
      pcri44[j]<-0
    }
  }
  
  pcri44 = data.frame(pcri44)
  prob4<-sum(pcri44==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+4]<DL
        & df_totalM$obs[i+8]<DL
        & df_totalM$obs[i+12]<DL
        #& df_totalM$obs[i+16]<DL 
        & df_totalM$time[i]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+8]
        & df_totalM$time[i+8]<df_totalM$time[i+12]
        #& df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      bcri4[i]<-df_totalM$time[i+12]
    }else {
      bcri4[i]<-30000
    }
  }
  
  bcri4 = data.frame(bcri4)
  df_totalM<-cbind(df_totalM,bcri4)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri44[j]<-min(df_totalM3$bcri4) - min(df_totalM3$cri)
    
  }
  
  bcri44 = data.frame(bcri44)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri44$pcri44[j]==1){
      lcri4[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri4)
    }else {
      lcri4[j]<-NaN
    }
  }
  
  lcri4 = data.frame(lcri4)
  lcri4$lcri4[lcri4$lcri4 < 0] <- 0
  lcri4 <- na.omit(lcri4)
  if (length(lcri4$lcri4) > 0) {
    lcri44 <- mean(lcri4$lcri4)
  }else {
    lcri44 <- 0
  }
  
  
  #################### Interval of tests = 5 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL 
        & df_totalM$obs[i+10]<DL 
        & df_totalM$obs[i+15]<DL 
        #& df_totalM$obs[i+20]<DL
        #& df_totalM$V[i+5]>DL
        #& df_totalM$V[i+10]>DL
        & df_totalM$V[i+15]>IF
        #& df_totalM$V[i+20]>DL
        & df_totalM$time[i]<df_totalM$time[i+5] 
        & df_totalM$time[i+5]<df_totalM$time[i+10] 
        & df_totalM$time[i+10]<df_totalM$time[i+15]
        #& df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      pcri5[i]<-1
    }else {
      pcri5[i]<-0
    }
  }
  
  pcri5 = data.frame(pcri5)
  df_totalM<-cbind(df_totalM,pcri5)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri5 >= 1) >=1){
      pcri55[j]<-1
    }else {
      pcri55[j]<-0
    }
  }
  
  pcri55 = data.frame(pcri55)
  prob5<-sum(pcri55==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL
        & df_totalM$obs[i+10]<DL
        & df_totalM$obs[i+15]<DL
        #& df_totalM$obs[i+20]<DL 
        & df_totalM$time[i]<df_totalM$time[i+5] 
        & df_totalM$time[i+5]<df_totalM$time[i+10]
        & df_totalM$time[i+10]<df_totalM$time[i+15]
        #& df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      bcri5[i]<-df_totalM$time[i+15]
    }else {
      bcri5[i]<-30000
    }
  }
  
  bcri5 = data.frame(bcri5)
  df_totalM<-cbind(df_totalM,bcri5)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri55[j]<-min(df_totalM3$bcri5) - min(df_totalM3$cri)
    
  }
  
  bcri55 = data.frame(bcri55)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri55$pcri55[j]==1){
      lcri5[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri5)
    }else {
      lcri5[j]<-NaN
    }
  }
  
  lcri5 = data.frame(lcri5)
  lcri5$lcri5[lcri5$lcri5 < 0] <- 0
  lcri5 <- na.omit(lcri5)
  if (length(lcri5$lcri5) > 0) {
    lcri55 <- mean(lcri5$lcri5)
  }else {
    lcri55 <- 0
  }
  
  
  #################### Result ####################
  Prob4[l,] <- c(prob1,prob2,prob3,prob4,prob5)
  
  B4 <- data.frame(bcri11,bcri22,bcri33,bcri44,bcri55)
  Burden4 <- rbind(Burden4,B4)
  
  L4 <- data.frame(lcri11,lcri22,lcri33,lcri44,lcri55)
  Length4 <- rbind(Length4,L4)
  
  
  df_totalM<-subset(df_totalM,select=-c(pcri1,pcri2,pcri3,pcri4,pcri5,bcri1,bcri2,bcri3,bcri4,bcri5))
  pcri1<-c(); pcri2<-c(); pcri3<-c(); pcri4<-c(); pcri5<-c();
  pcri11<-c(); pcri22<-c(); pcri33<-c(); pcri44<-c(); pcri55<-c();
  bcri1<-c(); bcri2<-c(); bcri3<-c(); bcri4<-c(); bcri5<-c();
  bcri11<-c(); bcri22<-c(); bcri33<-c(); bcri44<-c(); bcri55<-c();
  lcri1<-c(); lcri2<-c(); lcri3<-c(); lcri4<-c(); lcri5<-c();
  rm(df_totalM3)
  
  
  
  
  ###############################################################################################################
  ######################################## 5 Consecutive negative results #######################################
  ###############################################################################################################
  
  #################### Interval of tests = 1 day ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL 
        & df_totalM$obs[i+2]<DL 
        & df_totalM$obs[i+3]<DL 
        & df_totalM$obs[i+4]<DL
        #& df_totalM$V[i]>DL
        #& df_totalM$V[i+1]>DL 
        #& df_totalM$V[i+2]>DL 
        #& df_totalM$V[i+3]>DL
        & df_totalM$V[i+4]>IF 
        & df_totalM$time[i]<df_totalM$time[i+1] 
        & df_totalM$time[i+1]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+3]
        & df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      pcri1[i]<-1
    }else {
      pcri1[i]<-0
    }
  }
  
  pcri1 = data.frame(pcri1)
  df_totalM<-cbind(df_totalM,pcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri1 >= 1) >=1){
      pcri11[j]<-1
    }else {
      pcri11[j]<-0
    }
  }
  
  pcri11 = data.frame(pcri11)
  prob1<-sum(pcri11==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+1]<DL
        & df_totalM$obs[i+2]<DL
        & df_totalM$obs[i+3]<DL
        & df_totalM$obs[i+4]<DL
        & df_totalM$time[i]<df_totalM$time[i+1] 
        & df_totalM$time[i+1]<df_totalM$time[i+2]
        & df_totalM$time[i+2]<df_totalM$time[i+3]
        & df_totalM$time[i+3]<df_totalM$time[i+4]
    ){
      bcri1[i]<-df_totalM$time[i+4]
    }else {
      bcri1[i]<-30000
    }
  }
  
  bcri1 = data.frame(bcri1)
  df_totalM<-cbind(df_totalM,bcri1)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri11[j]<-min(df_totalM3$bcri1) - min(df_totalM3$cri)
    
  }
  
  bcri11 = data.frame(bcri11)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri11$pcri11[j]==1){
      lcri1[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri1)
    }else {
      lcri1[j]<-NaN
    }
  }
  
  lcri1 = data.frame(lcri1)
  lcri1$lcri1[lcri1$lcri1 < 0] <- 0
  lcri1 <- na.omit(lcri1)
  if (length(lcri1$lcri1) > 0) {
    lcri11 <- mean(lcri1$lcri1)
  }else {
    lcri11 <- 0
  }
  
  
  #################### Interval of tests = 2 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+2]<DL 
        & df_totalM$obs[i+4]<DL 
        & df_totalM$obs[i+6]<DL 
        & df_totalM$obs[i+8]<DL 
        #& df_totalM$V[i+2]>DL
        #& df_totalM$V[i+4]>DL
        #& df_totalM$V[i+6]>DL
        & df_totalM$V[i+8]>IF
        & df_totalM$time[i]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+6]
        & df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      pcri2[i]<-1
    }else {
      pcri2[i]<-0
    }
  }
  
  pcri2 = data.frame(pcri2)
  df_totalM<-cbind(df_totalM,pcri2)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri2 >= 1) >=1){
      pcri22[j]<-1
    }else {
      pcri22[j]<-0
    }
  }
  
  pcri22 = data.frame(pcri22)
  prob2<-sum(pcri22==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+2]<DL
        & df_totalM$obs[i+4]<DL
        & df_totalM$obs[i+6]<DL
        & df_totalM$obs[i+8]<DL
        & df_totalM$time[i]<df_totalM$time[i+2] 
        & df_totalM$time[i+2]<df_totalM$time[i+4]
        & df_totalM$time[i+4]<df_totalM$time[i+6]
        & df_totalM$time[i+6]<df_totalM$time[i+8]
    ){
      bcri2[i]<-df_totalM$time[i+8]
    }else {
      bcri2[i]<-30000
    }
  }
  
  bcri2 = data.frame(bcri2)
  df_totalM<-cbind(df_totalM,bcri2)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri22[j]<- min(df_totalM3$bcri2) - min(df_totalM3$cri)
    
  }
  
  bcri22 = data.frame(bcri22)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri22$pcri22[j]==1){
      lcri2[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri2)
    }else {
      lcri2[j]<-NaN
    }
  }
  
  lcri2 = data.frame(lcri2)
  lcri2$lcri2[lcri2$lcri2 < 0] <- 0
  lcri2 <- na.omit(lcri2)
  if (length(lcri2$lcri2) > 0) {
    lcri22 <- mean(lcri2$lcri2)
  }else {
    lcri22 <- 0
  }
  
  
  #################### Interval of tests = 3 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+3]<DL 
        & df_totalM$obs[i+6]<DL 
        & df_totalM$obs[i+9]<DL 
        & df_totalM$obs[i+12]<DL
        #& df_totalM$V[i+3]>DL
        #& df_totalM$V[i+6]>DL
        #& df_totalM$V[i+9]>DL
        & df_totalM$V[i+12]>IF
        & df_totalM$time[i]<df_totalM$time[i+3] 
        & df_totalM$time[i+3]<df_totalM$time[i+6] 
        & df_totalM$time[i+6]<df_totalM$time[i+9]
        & df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      pcri3[i]<-1
    }else {
      pcri3[i]<-0
    }
  }
  
  pcri3 = data.frame(pcri3)
  df_totalM<-cbind(df_totalM,pcri3)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri3 >= 1) >=1){
      pcri33[j]<-1
    }else {
      pcri33[j]<-0
    }
  }
  
  pcri33 = data.frame(pcri33)
  prob3<-sum(pcri33==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+3]<DL
        & df_totalM$obs[i+6]<DL
        & df_totalM$obs[i+9]<DL
        & df_totalM$obs[i+12]<DL
        & df_totalM$time[i]<df_totalM$time[i+3] 
        & df_totalM$time[i+3]<df_totalM$time[i+6]
        & df_totalM$time[i+6]<df_totalM$time[i+9]
        & df_totalM$time[i+9]<df_totalM$time[i+12]
    ){
      bcri3[i]<-df_totalM$time[i+12]
    }else {
      bcri3[i]<-30000
    }
  }
  
  bcri3 = data.frame(bcri3)
  df_totalM<-cbind(df_totalM,bcri3)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri33[j]<-min(df_totalM3$bcri3) - min(df_totalM3$cri)
    
  }
  
  bcri33 = data.frame(bcri33)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri33$pcri33[j]==1){
      lcri3[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri3)
    }else {
      lcri3[j]<-NaN
    }
  }
  
  lcri3 = data.frame(lcri3)
  lcri3$lcri3[lcri3$lcri3 < 0] <- 0
  lcri3 <- na.omit(lcri3)
  if (length(lcri3$lcri3) > 0) {
    lcri33 <- mean(lcri3$lcri3)
  }else {
    lcri33 <- 0
  }
  
  
  #################### Interval of tests = 4 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL 
        & df_totalM$obs[i+4]<DL 
        & df_totalM$obs[i+8]<DL 
        & df_totalM$obs[i+12]<DL 
        & df_totalM$obs[i+16]<DL
        #& df_totalM$V[i+4]>DL
        #& df_totalM$V[i+8]>DL
        #& df_totalM$V[i+12]>DL
        & df_totalM$V[i+16]>IF
        & df_totalM$time[i]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+8] 
        & df_totalM$time[i+8]<df_totalM$time[i+12]
        & df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      pcri4[i]<-1
    }else {
      pcri4[i]<-0
    }
  }
  
  pcri4 = data.frame(pcri4)
  df_totalM<-cbind(df_totalM,pcri4)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri4 >= 1) >=1){
      pcri44[j]<-1
    }else {
      pcri44[j]<-0
    }
  }
  
  pcri44 = data.frame(pcri44)
  prob4<-sum(pcri44==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+4]<DL
        & df_totalM$obs[i+8]<DL
        & df_totalM$obs[i+12]<DL
        & df_totalM$obs[i+16]<DL
        & df_totalM$time[i]<df_totalM$time[i+4] 
        & df_totalM$time[i+4]<df_totalM$time[i+8]
        & df_totalM$time[i+8]<df_totalM$time[i+12]
        & df_totalM$time[i+12]<df_totalM$time[i+16]
    ){
      bcri4[i]<-df_totalM$time[i+16]
    }else {
      bcri4[i]<-30000
    }
  }
  
  bcri4 = data.frame(bcri4)
  df_totalM<-cbind(df_totalM,bcri4)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri44[j]<-min(df_totalM3$bcri4) - min(df_totalM3$cri)
    
  }
  
  bcri44 = data.frame(bcri44)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri44$pcri44[j]==1){
      lcri4[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri4)
    }else {
      lcri4[j]<-NaN
    }
  }
  
  lcri4 = data.frame(lcri4)
  lcri4$lcri4[lcri4$lcri4 < 0] <- 0
  lcri4 <- na.omit(lcri4)
  if (length(lcri4$lcri4) > 0) {
    lcri44 <- mean(lcri4$lcri4)
  }else {
    lcri44 <- 0
  }
  
  
  #################### Interval of tests = 5 days ####################
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL 
        & df_totalM$obs[i+10]<DL 
        & df_totalM$obs[i+15]<DL 
        & df_totalM$obs[i+20]<DL
        #& df_totalM$V[i+5]>DL
        #& df_totalM$V[i+10]>DL
        #& df_totalM$V[i+15]>DL
        & df_totalM$V[i+20]>IF
        & df_totalM$time[i]<df_totalM$time[i+5] 
        & df_totalM$time[i+5]<df_totalM$time[i+10] 
        & df_totalM$time[i+10]<df_totalM$time[i+15]
        & df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      pcri5[i]<-1
    }else {
      pcri5[i]<-0
    }
  }
  
  pcri5 = data.frame(pcri5)
  df_totalM<-cbind(df_totalM,pcri5)
  
  for( j in 1:N ) {
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (sum(df_totalM3$pcri5 >= 1) >=1){
      pcri55[j]<-1
    }else {
      pcri55[j]<-0
    }
  }
  
  pcri55 = data.frame(pcri55)
  prob5<-sum(pcri55==1)/N
  
  for( i in 1:TM ) {
    if (df_totalM$time[i]>=First
        & df_totalM$obs[i]<DL
        & df_totalM$obs[i+5]<DL
        & df_totalM$obs[i+10]<DL
        & df_totalM$obs[i+15]<DL
        & df_totalM$obs[i+20]<DL
        & df_totalM$time[i]<df_totalM$time[i+5] 
        & df_totalM$time[i+5]<df_totalM$time[i+10]
        & df_totalM$time[i+10]<df_totalM$time[i+15]
        & df_totalM$time[i+15]<df_totalM$time[i+20]
    ){
      bcri5[i]<-df_totalM$time[i+20]
    }else {
      bcri5[i]<-30000
    }
  }
  
  bcri5 = data.frame(bcri5)
  df_totalM<-cbind(df_totalM,bcri5)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    bcri55[j]<-min(df_totalM3$bcri5) - min(df_totalM3$cri)
    
  }
  
  bcri55 = data.frame(bcri55)
  
  for( j in 1:N ) {
    
    df_totalM3<-subset(df_totalM,df_totalM$ID==j)
    
    if (pcri55$pcri55[j]==1){
      lcri5[j]<-min(df_totalM3$cri) - min(df_totalM3$bcri5)
    }else {
      lcri5[j]<-NaN
    }
  }
  
  lcri5 = data.frame(lcri5)
  lcri5$lcri5[lcri5$lcri5 < 0] <- 0
  lcri5 <- na.omit(lcri5)
  if (length(lcri5$lcri5) > 0) {
    lcri55 <- mean(lcri5$lcri5)
  }else {
    lcri55 <- 0
  }
  
  
  #################### Result ####################
  Prob5[l,] <- c(prob1,prob2,prob3,prob4,prob5)
  
  B5 <- data.frame(bcri11,bcri22,bcri33,bcri44,bcri55)
  Burden5 <- rbind(Burden5,B5)
  
  L5 <- data.frame(lcri11,lcri22,lcri33,lcri44,lcri55)
  Length5 <- rbind(Length5,L5)
  
  
  df_totalM<-subset(df_totalM,select=-c(pcri1,pcri2,pcri3,pcri4,pcri5,bcri1,bcri2,bcri3,bcri4,bcri5))
  pcri1<-c(); pcri2<-c(); pcri3<-c(); pcri4<-c(); pcri5<-c();
  pcri11<-c(); pcri22<-c(); pcri33<-c(); pcri44<-c(); pcri55<-c();
  bcri1<-c(); bcri2<-c(); bcri3<-c(); bcri4<-c(); bcri5<-c();
  bcri11<-c(); bcri22<-c(); bcri33<-c(); bcri44<-c(); bcri55<-c();
  lcri1<-c(); lcri2<-c(); lcri3<-c(); lcri4<-c(); lcri5<-c();
  rm(df_totalM3)
  
  
}

Probability <- data.frame(Prob1,Prob2,Prob3,Prob4,Prob5)
colnames(Probability) <- c("p11","p12","p13","p14","p15",
                           "p21","p22","p23","p24","p25",
                           "p31","p32","p33","p34","p35",
                           "p41","p42","p43","p44","p45",
                           "p51","p52","p53","p54","p55")

Length <- data.frame(Length1,Length2,Length3,Length4,Length5)
colnames(Length) <- c("l11","l12","l13","l14","l15",
                      "l21","l22","l23","l24","l25",
                      "l31","l32","l33","l34","l35",
                      "l41","l42","l43","l44","l45",
                      "l51","l52","l53","l54","l55")

Burden <- data.frame(Burden1,Burden2,Burden3,Burden4,Burden5)
colnames(Burden) <- c("b11","b12","b13","b14","b15",
                      "b21","b22","b23","b24","b25",
                      "b31","b32","b33","b34","b35",
                      "b41","b42","b43","b44","b45",
                      "b51","b52","b53","b54","b55")

write_xlsx(Probability,"Probability_S61.xlsx")
write_xlsx(Length,"Length_S61.xlsx")
write_xlsx(Burden,"Burden_S61.xlsx")




