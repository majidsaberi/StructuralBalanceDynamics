rm(list = ls())
library(parallel)
setwd("/Users/majid/Projects/Frustration-Lifetime/Analysis/")  

load("Data-50.RData")

Pheno_unrel <- read.csv("HCP-Pheno-Unrelated.csv")
dim(TS_50)
TS_50 <- TS_50[match(Pheno_unrel$Subject,Pheno$Subject),,]

DynConn <- function(TS,WL=50){
  dm <- dim(TS)
  dconn <- array(dim = c(dm[2]-WL,dm[1],dm[1]) )
  for(i in 1:(dm[2]-WL)){
    dconn[i,,] <- cor(t(TS[,i:(i+WL)]),use = "na.or.complete")
  }
  return(dconn)
}

DynTriad <- function(dconn){
  dm <- dim(dconn)
  dconn[dconn < 0] <- -1
  dconn[dconn > 0] <- 1
  dtriad <- array(dim = c(dm[1],dm[2],dm[2],dm[2]))
  for(t in 1:dm[1]){
    for(i in 1:(dm[2]-2) ){
      for(j in (i+1):(dm[2]-1)){
        for(k in (j+1):(dm[2])){
          dtriad[t,i,j,k] <- dconn[t,i,j]+dconn[t,i,k]+dconn[t,j,k]
        }
      }
    }
  }
  return(dtriad)
}


LifeTime <- function(dtriad){
  dm <- dim(dtriad)
  lf <- array(dim = c(dm[2],dm[2],dm[2]) ) 
  lf <- array(list(),dim = c(dm[2],dm[2],dm[2]) ) 
  for(i in 1:(dm[2]-2)){
    for(j in (i+1):(dm[2]-1)){
      for(k in (j+1):(dm[2])){
        trts <- dtriad[,i,j,k]
        lf[[i,j,k]] <- rle(trts)
      }
    }
  }
  return(lf)
}


LF_WholeBrain <- function(lf){
  neg1 <- apply(lf,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == -1],na.rm=T) )
  pos3 <- apply(lf,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == 3],na.rm=T) )
  pos1 <- apply(lf,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == 1],na.rm=T) )
  neg3 <- apply(lf,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == -3],na.rm=T) )
  
  output <- c(mean(neg1,na.rm=T),
              mean(pos3,na.rm=T),
              mean(pos1,na.rm=T),
              mean(neg3,na.rm=T) )
  
  names(output) <- c("Neg1","Pos3","Pos1","Neg3")
  
  return(output)
}

DynTriadEnergy <- function(dconn){
  dm <- dim(dconn)
  dtriadenergy <- array(dim = c(dm[1],dm[2],dm[2],dm[2]))
  for(t in 1:dm[1]){
    for(i in 1:(dm[2]-2) ){
      for(j in (i+1):(dm[2]-1)){
        for(k in (j+1):(dm[2])){
          a <- - dconn[t,i,j]*dconn[t,i,k]*dconn[t,j,k]
          dtriadenergy[t,i,j,k]  <- sign(a) * abs(a)^(1/3)
        }
      }
    }
  }
  return(dtriadenergy)
}

DynTriadPeakEnergy <- function(dtriadenergy,lf){
  dm <- dim(dtriadenergy)
  de <- lf
  for(i in 1:(dm[2]-2)){
    for(j in (i+1):(dm[2]-1)){
      for(k in (j+1):(dm[2])){
        trts1 <- dtriadenergy[,i,j,k]
        trts2 <- lf[[i,j,k]]
        kk <- 1
        for(jj in 1:length(trts2$lengths)){
          a <- trts1[kk:(kk+(trts2$lengths)[jj]-1)]
          de[[i,j,k]]$lengths[jj] <- max(abs(a))
          kk <- kk + (trts2$lengths)[jj] 
        }
      }
    }
  }
  return(de)
}        

DPE_WholeBrain <- function(de){
  neg1 <- apply(de,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == -1],na.rm=T) )
  pos3 <- apply(de,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == 3],na.rm=T) )
  pos1 <- apply(de,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == 1],na.rm=T) )
  neg3 <- apply(de,c(1,2,3),function(x) mean(x[[1]][[1]][ x[[1]][[2]] == -3],na.rm=T) )
  
  output <- c(mean(neg1,na.rm=T),
              mean(pos3,na.rm=T),
              mean(pos1,na.rm=T),
              mean(neg3,na.rm=T) )
  
  names(output) <- c("Neg1","Pos3","Pos1","Neg3")
  
  return(output)
}

#


###
##
#

mainfunction <- function(subjnum){
  dconn <- DynConn(TS_50[subjnum,,],WL = 50)
  
  dtriad <- DynTriad(dconn)
  lf <- LifeTime(dtriad) 
  lfb <- LF_WholeBrain(lf)
  
  dtriadenergy <- DynTriadEnergy(dconn)
  de <- DynTriadPeakEnergy(dtriadenergy,lf)
  deb <- DPE_WholeBrain(de)  

  trb <- c( sum(dtriad == -1 ,na.rm = T),
     sum(dtriad == 3 ,na.rm = T),
     sum(dtriad == 1 ,na.rm = T),
     sum(dtriad == -3 ,na.rm = T)) / sum(!is.na(dtriad))
  names(trb) <- names(lfb)
  
  output <- list(lfb,deb,trb)
  names(output) <- c("lifetime","peakenergy","traidnum")
  return(output)
}

##
#

alllifetime <- list()

for(i in 1:100){
  alllifetime[[i]] <- mainfunction(i)
  gc()
  print(i)
}

alllifetime[[100]]

rm(list = ls()[!ls() %in% c("alllifetime","Pheno_unrel")])

