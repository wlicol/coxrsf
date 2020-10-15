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
# 1. D-index/RSF (OOB)
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
# rfsrc.obj$yvar is the dependent variable of the form Surv(Time, Status) 
#
# If ave.flag=1 (default), the output is the averaged ("I" in IBS) Brier score, without any weight
# the survival function is only available from rfsrc.obj at "unique dead times".
#
# If ave.flag=0, a matrix with three rows are the output. The first row is IBS for people who
# are dead (status=1); the second row is IBS for people who are censored (status=0); and the
# last row the average of the two groups.
#
# Note that our calculation relies on object from rfsrc, without getting into each bootstrap+OOB split
# of samples in one tree, and we trust each sample should have a chance to be in OOB, and
# rfsrc.obj$survival.oob provides an average of all trees.
############################################################################################

ibsRSFoob <- function(rfsrc.obj, ave.flag=1){
 nsample <- rfsrc.obj$n
 time <- rfsrc.obj$time.interest
 ntime <- length(time)

 deadStatus <- rfsrc.obj$yvar[,2]
 deadtime <- rfsrc.obj$yvar[deadStatus==1,1]	# may not be identical to sort(time), but = unique( sort(time))
 censortime <- rfsrc.obj$yvar[deadStatus==0,1]

 survP <- rfsrc.obj$survival.oob
 if( length(survP) == (nsample*ntime)){
  tmp <- matrix(survP, ncol=ntime)
  tmp1 <- tmp[deadStatus==1, ]
  tmp0 <- tmp[deadStatus==0, ]
  bs1 <- matrix(nrow=nrow(tmp1), ncol=ntime)
  bs0 <- matrix(nrow=nrow(tmp0), ncol=ntime)

  for(idead in 1:nrow(tmp1)){
   coreBefore <- c(1:ntime)[ time < deadTime[idead] ]
   coreAfter <- c(1:ntime)[ time >= deadTime[idead] ]
   bs1[idead, coreBefore] <- (1-tmp1[idead,coreBefore])^2
   bs1[idead, coreAfter] <- tmp1[idead,coreAfter]^2
  }

  for(ilive in 1:nrow(tmp0)){
   coreBefore <- c(1:ntime)[ time < censorTime[ilive] ]
   coreAfter <- c(1:ntime)[ time >= censorTime[ilive] ]
   bs0[ilive, coreBefore] <- (1-tmp0[ilive,coreBefore])^2
   bs0[ilive, coreAfter] <-  NA			# it should not be zero.
  }

  bs <- rbind(bs1, bs0)
  ibs.time1 <- apply(bs1, 2, mean)
  ibs.time0 <- apply(bs0, 2, function(x){ mean( sort(x))  })
  ibs.time <- apply(bs, 2, function(x){ mean( sort(x))  } )
  ibs.time.all <- rbind(ibs.time1,ibs.time0,ibs.time)

  if(ave.flag==1){
	aveibs <- mean(ibs.time.all[3,])
	return( aveibs)
  } else if(ave.flag==0){
	return( ibs.time.all)
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
############################################################################################

library(survival)

dindexCox <- function(cox.obj, newx, newy){

 pred.risk.x  <- as.numeric(predict(cox.obj, newdata=newx, type="risk"))
 time.y <- newy[, 1]
 status.y <- newy[, 2]
 nsample <- length(time.y)

 ord <- order(time.y, -status.y)
 time <- time.y[ord]		# this still may be identical to unique(time)
 status <- status.y[ord]
 pred.risk <- pred.risk.x[ord]

 tmp0 <- cbind(1:nsample, pred.risk, time,status)
 tmp <- tmp0[sort.list(tmp0[,3]), ]

 wh <- c(1:nsample)[tmp[,4]==1]
 total <- 0
 concordance <- 0

 # the first person i is the dead person
 for(i in wh){
  if(i < nsample){
   # the second person j can be either dead or cencored
   for(j in (i+1):nsample){

    #person-j live longer
    if(tmp[j,3] > tmp[i,3]){
     total <- total+1
     #person-j also has a lower (predicted) risk, so concordance 
     if(tmp[j,2] < tmp[i,2]){
	concordance <- concordance +1
     # in case the risks are tied
     } else if(tmp[j,2]==tmp[i,2]){
	concordance <- concordance +0.5
     }
  }
 }
 }
 }
 return( list( concordant=concordant, total=total, cindex=concordant/total, dindex= 1-concordant/total ))


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

library(survival)


ibsCox2 <- function(cox.obj, newx, newy, ave.flag=1){
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
 bs1 <- matrix(nrow=nrow(tmp1), ncol=ntime)
 bs0 <- matrix(nrow=nrow(tmp0), ncol=ntime)

 for(idead in 0:nrow(tmp1)){
  coreBefore <- c(1:ntime)[ time < deadTime[idead] ]
  coreAfter <- c(1:ntime)[ time >= deadTime[idead] ]
  bs1[idead, coreBefore] <- (1-tmp1[idead,coreBefore])^2
  bs1[idead, coreAfter] <- tmp1[idead,coreAfter]^2
 }

 for(ilive in 1:nrow(tmp0)){
  coreBefore <- c(1:ntime)[ time < censorTime[ilive] ]
  coreAfter <- c(1:ntime)[ time >= censorTime[ilive] ]
  bs0[ilive, coreBefore] <- (1-tmp0[ilive,coreBefore])^2
  bs0[ilive, coreAfter] <-  NA			# it should not be zero.
 }

 bs <- rbind(bs1, bs0)
 ibs.time1 <- apply(bs1, 2, mean)
 ibs.time0 <- apply(bs0, 2, function(x){ mean( sort(x))  })
 ibs.time <- apply(bs, 2, function(x){ mean( sort(x))  } )
 ibs.time.all <- rbind(ibs.time1,ibs.time0,ibs.time)

 if(ave.flag==1){
	aveibs <- mean(ibs.time.all[3,])
	return( aveibs)
 } else if(ave.flag==0){
	return( ibs.time.all)
 }
}

