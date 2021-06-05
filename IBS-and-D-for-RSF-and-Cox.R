############################################################################################
# We provide simple R code to calculate OOB (out of bag)
# errors, either D-index =1-Cindex or IBS (integrated Brier's score),
# for both RSF (random survival forest) object generated
# by rfsrc() from randomForestSRC package and by coxph() from
# survival package.
#
# randomForestSRC: https://cran.r-project.org/web/packages/randomForestSRC/index.html
# survival: https://cran.r-project.org/web/packages/survival/index.html
#
# we call the object created by rfsrc "rfsrc.obj": rfsrc.obj <- rfsrc(....)
# we call the object created by coxph "cox.obj":  cox.obj    <- coxph(....)
#
# Since the "C" in C-index means concordance, we call 1-Cindex "D-index" for
# discordance index. (Harrell (1982), Evaluating the yield of medical tests, 247:2543-2546.)
#
# The IBS (integrated Brier's score is an average of Brier's score (Brier (1950), Verification 
# of forecasts expressed in terms of probability,, Monthly Weather Rev. 78:1-3) at time points,
# unweighted in our implementation.
#
############################################################################################


############################################################################################
# 1. D-index/RSF
#
# This error for OOB is already provided by rfsrc.obj:
#	rfsrc.obj$err.rate[rfsrc.obj$ntree]
# See this description from rfsrc manual page:
#  "Prediction error is calculated using OOB data...For survival, prediction error is measured 
# by 1-C, where C is Harrell’s (Harrell et al., 1982) concordance index."
#
# For a typical run, rfsrc.obj$err.rate returns an array which are all NA except the last element.
############################################################################################


############################################################################################
# 2. IBS/RSF (OOB) 
#
# rfsrc.obj$time.interest lists the "ordered unique death time"
# rfsrc.obj$survival.oob is the "OOB survival function", for all samples, at all time of death,
# as an array. The order of these suvival functions are first person's survival function at
# the sorted time point, followed by the second person's, etc.
#
# if self=1, the "in-the-bag" obj$survival is used (default value is self=0, or use OOB samples).
#
# rfsrc.obj$yvar is the dependent variable of the form Surv(Time, Status) 
#
# If ave.flag=0,  the function returns a three-row, T-column (T is the number of unique death time) matrix.
# The first row is the Brier Score for status=1 samples, at each time point.
# The second row is the Brier Score for status=0 samples, at each time point.
# The third row is the Brier Score for all samples, at each time point.
#
# If ave.flag=1 (default), we have the integrated BS (IBS), i.e., mean(ibs.time.all[3,])
#
# Note that our calculation relies on object from rfsrc, without getting into each bootstrap+OOB split
# of samples in one tree, and we trust each sample should have a chance to be in OOB, and
# rfsrc.obj$survival.oob provides an average of all trees.
############################################################################################


ibsRSF <- function(rfsrc.obj, self=0, ave.flag=1){
 obj <- rfsrc.obj
 nsample <- obj$n
 time <- obj$time.interest
 ntime <- length(time)

 deadStatus <- obj$yvar[,2]
 deadTime <- obj$yvar[deadStatus==1,1]
 censorTime <- obj$yvar[deadStatus==0,1]

 if(self==0){
  survP <- obj$survival.oob
 }
 else{
  if(self==1){
   survP <- obj$survival
 }}

if(length(survP) == (nsample*ntime)){

   tmp <- matrix(survP, ncol=ntime)
   tmp1 <- tmp[deadStatus==1, ]
   tmp0 <- tmp[deadStatus==0, ]
   bs1 <- matrix(nrow=nrow(tmp1), ncol=ntime)
   bs0 <- matrix(nrow=nrow(tmp0), ncol=ntime)

   # go through dead 
   for(si in 1:nrow(tmp1)){
     coreBefore <- c(1:ntime)[ time < deadTime[si] ]
     coreAfter <- c(1:ntime)[ time >= deadTime[si] ]
     bs1[si, coreBefore] <- (1-tmp1[si,coreBefore])^2
     bs1[si, coreAfter] <- tmp1[si,coreAfter]^2
   }

   # go through censored
   for(si in 1:nrow(tmp0)){
     coreBefore <- c(1:ntime)[ time < censorTime[si] ]
     coreAfter <- c(1:ntime)[ time >= censorTime[si] ]
     bs0[si, coreBefore] <- (1-tmp0[si,coreBefore])^2
     # it should not be zero, but unknown (not contributing)
     # bs0[si, coreAfter] <- 0
     bs0[si, coreAfter] <- NA
   }

  bs <- rbind(bs1, bs0)
  ibs.time1 <- apply(bs1, 2, mean)
  ibs.time0 <- apply(bs0, 2, function(x){ mean( sort(x))  })
  ibs.time <- apply(bs, 2, function(x){ mean( sort(x))  } )
  ibs.time.all <- rbind(ibs.time1,ibs.time0,ibs.time)

  if(ave.flag==1){
   aveibs <- mean(ibs.time.all[3,])
   return(aveibs)
  } else if(ave.flag==0){
   return(ibs.time.all)
  }

}}



############################################################################################
# 3. D-index/Cox on new data
# 
# Unlike IBS, there is no survival function at time points, only a single predicted risk per person.
#
# To get OOB D-index, the bootstrap+OOB split is done outside this function. Then the function is
# applied to the OOB set.
#
# predict from "survival" package is needed
#
# New data is split between independent variables (newx) and dependent variable (newy) which
# of the form Surv(time, status). newx is used to predict, newy is used to evaluate errors.
#
# A rarely used option is self=1, then newx/newy is not used. The IBS is calculated on the training dataset.
#
# The output contains these information: (1) total number of comparison; (2) number of concordance
# out of these comparisons; (3) C-index = #(concordance)/#(comparison); (4) D-index= 1-Cindex
############################################################################################

# library(survival)

dindexCOX <- function(cox.obj, newx, newy, self=0){

 # self=0: use newx to predict the risk, newy to calculate c-index
 if(self == 0){
  # default is log(risk)
  pred.risk  <- as.numeric(predict(cox.obj, newdata=newx, type="risk"))
  p <- ncol(newy)
  time <- newy[ ,p-1]
  status <- newy[ ,p]
 # self=1:  risk on the training set itself, original y (extracted from cox.obj)  to calculate c-index
 } else{
  pred.risk <- as.numeric(predict(cox.obj, type="risk"))
  oriy <- cox.obj$y
  p <- ncol(oriy)
  time <- oriy[ ,p-1]
  status <- oriy[ ,p]
 }

 n <- length(time)
 ord <- order(time, -status)
 time.ordered <-  time[ord]
 status.ordered <-  status[ord]
 pred.risk.ordered <- pred.risk[ord]

 tmp0 <- cbind(1:n, pred.risk, time,status)
 tmp.order <- tmp0[sort.list(tmp0[, 3]), ]

 wh <- c(1:n)[tmp.order[,4]==1]
 # wh <- which(status == 1)
 total <- 0
 concordant <- 0

 for(i in wh){
  # last point will not form a pair
  if(i < n){
   for(j in (i+1):n ){
     if(tmp.order[j,3] > tmp.order[i,3]){
        # in case of bootstrap, avoid compare sample with itself
        total <- total+1
        if(tmp.order[j,2] < tmp.order[i,2]){
         concordant <- concordant +1
	}
	else{ 
         if(tmp.order[j,2]==tmp.order[i,2]){
          concordant <- concordant + 0.5
	 }
	}
     }
   } # end of j loop
  } # end of condition i < n (or i!=n)
 } # end of i loop
 return( list(concordant=concordant, total=total, cindex=concordant/total, dindex= 1-concordant/total ))

}



############################################################################################
# 4. IBS/Cox on new data
#
# surfit from survival R package is needed
#
# The predicted survival function at time points are:
#  survfit(cox.obj, newx)$surv 
# for "the estimate of survival at time t+0"
############################################################################################

# library(survival)


ibsCOX  <- function(cox.obj, newx, newy, ave.flag=1){
 # (new)data information is derived from newy, not from cox.obj, which is derived from a different (part of the) dataset
 deadStatus <- newy[,2]
 deadTime <- newy[deadStatus==1,1]
 censorTime <- newy[deadStatus==0,1]


 survfit.result <- survfit(cox.obj, newdata=newx)
 # NOTE: the time from rfsrc.obj$time.interest only contains "dead time". here it contains dead/censored times
 time <- survfit.result$time
 ntime <- length(time)
 # survival function matrix from survfit() has ntime rows, nsample columns. to be consistent with ibsRSFoob, it's transposed
 tmp <- t( survfit.result$surv )

 # the remaining part is identical to ibsRSFoob
 tmp1 <- tmp[deadStatus==1, ]
 tmp0 <- tmp[deadStatus==0, ]

# what if there is only one status=1 sample
if( length(deadTime) > 1){
 bs1 <- matrix(nrow=nrow(tmp1), ncol=ntime)

 for(idead in 0:nrow(tmp1)){
  coreBefore <- c(1:ntime)[ time < deadTime[idead] ]
  coreAfter <- c(1:ntime)[ time >= deadTime[idead] ]
  bs1[idead, coreBefore] <- (1-tmp1[idead,coreBefore])^2
  bs1[idead, coreAfter] <- tmp1[idead,coreAfter]^2
 }
}

# what if there is only one status=0 sample
if( length(censorTime) > 1){
 bs0 <- matrix(nrow=nrow(tmp0), ncol=ntime)

 for(ilive in 1:nrow(tmp0)){
  coreBefore <- c(1:ntime)[ time < censorTime[ilive] ]
  coreAfter <- c(1:ntime)[ time >= censorTime[ilive] ]
  bs0[ilive, coreBefore] <- (1-tmp0[ilive,coreBefore])^2
  bs0[ilive, coreAfter] <-  NA			# it should not be zero.
 }
}

# what if there is only one status=0 or one status=1 sample
if(length(deadTime) > 1 & length(censorTime) > 1 ){

 bs <- rbind(bs1, bs0)
 ibs.time1 <- apply(bs1, 2, mean)
 ibs.time0 <- apply(bs0, 2, function(x){ mean( sort(x))  })
 ibs.time <- apply(bs, 2, function(x){ mean( sort(x))  } )
 ibs.time.all <- rbind(ibs.time1,ibs.time0,ibs.time)
}else{
 if(length(censorTime) > 1){
  ibs.time.all <- rbind( rep(NA, length(ibs.time0)), ibs.time0, ibs.time0)
 }
 if(length(deadTime) > 1){
  ibs.time.all <- rbind(ibs.time1, rep(NA, length(ibs.time1)), ibs.time1)
 }
}

 if(ave.flag==1){
	aveibs <- mean(ibs.time.all[3,])
	return( aveibs)
 } else if(ave.flag==0){
	return( ibs.time.all)
 }
}

