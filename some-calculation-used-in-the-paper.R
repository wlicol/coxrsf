#-----------------------------------------------------------
# please consult the paper:
#
# Sirin Cetin, Ayse Ulgen, Isa Dede, Wentian Li (2021),
# "On fair performance comparison between random survival forest 
# and Cox regression: an example of colorectal cancer study"
# SciMedicine Journal, 3(1):66-76.
#
# doi: 10.28991/SciMedJ-2021-0301-9
#-----------------------------------------------------------

cr <- read.table("colorectal.txt", header=T)
names(cr) <- c("age", "gender", "ctype", "leukocyte", "neutrophile","lymphocyte", "platelet","status", "time")



library(randomForestSRC)
# use the RSF function rfsrc
library(survival)
# need the Surv function

source("IBS-and-D-for-RSF-and-Cox.R")
# use the Dindex and IBS error calculation

#####################################
# calculations for Fig.1: Dindex or IBS 
# as a function of ntree and nodesize
#####################################

# nrep <- 1000
nrep <- 10
res.par <- matrix(nrow=nrep, ncol=4)

for(i in 1:nrep){

 # ntree: range from 2-2000. default=1000
 ntree_tmp <- sample(c(2:2000))[1] #default: 1000

 # nodesize: range from 1-20. default=15
 nodesize_tmp <- sample(1:20)[1] 

 rfsrc.obj <- rfsrc(Surv(time,status) ~ ., data=cr, nodesize=nodesize_tmp, ntree=ntree_tmp, samptype="swr")

 # from rfsrc object directly.  Dindex=1-Cindex
 res.par[i,1]  <- rfsrc.obj$err[ length(rfsrc.obj$err)]
 # our own IBS function
 res.par[i,2]  <- ibsRSF(rfsrc.obj)

 res.par[i,3] <- ntree_tmp
 res.par[i,4] <- nodesize_tmp

}

# res.par matrix contains the information needed for Fig.1

#########################################
# the "0.632" IBS error calculation for 
# Fig.2  (for RSF)
#########################################

# weight for out-of-bag (=0.632)
wOOB <- 1-exp(-1)
# weight for in-the-bag (=0.368)
wIB <- exp(-1)

IBS.RSF.632 <- NULL
nite <- 100
for(i in 1:nite){
 rfsrc.obj <- rfsrc(Surv(time, status) ~ ., data=cr, ntree=500, samptype = "swr")
 err1.tmp <- ibsRSF(rfsrc.obj)
 err0.tmp <- ibsRSF(rfsrc.obj,self=1)
 errA.tmp <- wOOB*err1.tmp + wIB*err0.tmp
 IBS.RSF.632 <- c(IBS.RSF.632, errA.tmp)
}

#########################################
# the "0.632" IBS error calculation for 
# Fig.2  (for COX)
#########################################

# weight for out-of-bag (=0.632)
wOOB <- 1-exp(-1)
# weight for in-the-bag (=0.368)
wIB <- exp(-1)

IBS.COX.632 <- NULL
nite <- 100
for(i in 1:nite){

  ibs_tmp1 <- NULL
  ibs_tmp0 <- NULL

  # the number 500 is similar to ntree in RSF
  # for(j in 1:500){
  for(j in 1:5){

   train <- sample(1:nrow(cr), nrow(cr), replace=TRUE)

   cr.train <- cr[train, ] 
   cox.obj <- coxph(Surv(time, status)~ ., cr.train)
   newy0 <- Surv(cr.train$time, cr.train$status)
   newx0 <- subset(cr.train, select=-c(time, status))
   err.tmp0 <- ibsCOX(cox.obj, newx0, newy0)
   ibs_tmp0 <- c(ibs_tmp0, err.tmp0)

   cr.val <- cr[-train, ] 
   newy1 <- Surv(cr.val$time, cr.val$status)
   newx1 <- subset(cr.val, select=-c(time, status))
   err.tmp1 <- ibsCOX(cox.obj, newx1, newy1)
   ibs_tmp1 <- c(ibs_tmp1, err.tmp1)

  }
 IBS.COX.632 <- c(IBS.COX.632, mean(wOOB*ibs_tmp1+wIB*ibs_tmp0)) 
}


#########################################
# the 10-fold IBS error calculation for 
# Fig.2  (for COX)
#########################################

k <- 10
lastP <- nrow(cr)
border <- seq(1, lastP, as.integer(lastP/k))

binStart <- border[ -length(border)] 
binEnd <-  border[-1]-1
binEnd[ length(binEnd)] <- lastP

IBS.COX.10f <- NULL
# nite <- 100
nite <- 5
for(i in 1:nite){

 randomi <- sample(1:lastP)

 err_tmp <- NULL
 for(j in 1:k){
   validS <- randomi[binStart[j]:binEnd[j]]
   trainS <- c(1:lastP)[ -validS]
   cox.obj <- coxph(Surv(time, status) ~ ., cr[trainS,])

   newx <- subset(cr[validS,], select=-c(time, status) )
   newy <- Surv(cr[validS,]$time, cr[validS,]$status)
   err_tmp <- c(err_tmp, ibsCOX(cox.obj, newx, newy) ) 
 }
 # mean of the 10  10-fold errors
 IBS.COX.10f <- c(IBS.COX.10f, mean(err_tmp))
}




