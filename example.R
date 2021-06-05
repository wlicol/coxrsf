#------------------------------------------------------------------------------------------------------------
# this is a R script to illustrate the use of 
# (1) ibsRSF (rfsrc.obj, self=0,ave.flag=1) : use rfsrc.obj to calculate IBS error for out-of-bag (OOB) samples
# (2) dindexCOX(cox.obj, newx, newy,self=0) : use cox.obj to calculate Dindex error for new samples
# (3) ibsCOX( cox.obj, newx, newy, ave.flag=1) : use cox.obj to calculate IBS error for new samples
#
# the dataset used is 221-person colorectal cancer survival data ("colorectal.txt")
#
# Reference:
#
# Sirin Cetin, Ayse Ulgen, Isa Dede, W Li (2021),
# On fair performance comparison between random survival forest and Cox regression: an example of colorectal cancer study 
# SciMedicine Journal, 3(1):66-67.
#------------------------------------------------------------------------------------------------------------


cr <- read.table("colorectal.txt", header=T)
# make column header name simpler
names(cr) <- c("age", "gender", "ctype", "leukocyte", "neutrophile","lymphocyte", "platelet","status", "time")


library(survival)
library(randomForestSRC)

# our three functions are in this file
source("IBS-and-D-for-RSF-and-Cox.R")

surv.f  <- as.formula( Surv(time, status) ~ .)
Ntree <- 500
Nnode <- 15

#---------------------------------------------------------------------
# task-1. random survival forest, get D-index error.  do it 100 times
# default sample is "sampling without replacement" (swor).  to be 
# consistent with the original RF, use "sampling with replacement" (swr)
#---------------------------------------------------------------------

D.RSF.oob <- NULL
for(i in 1:100){
 rfsrc.obj <- rfsrc(surv.f, data=cr, nodesize=Nnode, ntree=Ntree, samptype="swr")
 D.RSF.oob <- c(D.RSF.oob, rfsrc.obj$err.rate[ rfsrc.obj$ntree] )
}

# one typical run
# sort(D.RSF.oob)
#  [1] 0.4395194 0.4396994 0.4398794 0.4399694 0.4401494 0.4407794 0.4411394
#  [8] 0.4415894 0.4419494 0.4421294 0.4426694 0.4428494 0.4433894 0.4435694
# [15] 0.4438394 0.4440194 0.4445594 0.4448294 0.4449194 0.4449194 0.4451895
# [22] 0.4452795 0.4453695 0.4454595 0.4456395 0.4457295 0.4458195 0.4459995
# [29] 0.4459995 0.4459995 0.4461795 0.4469895 0.4472595 0.4474395 0.4475295
# [36] 0.4475295 0.4475295 0.4476195 0.4476195 0.4477095 0.4477095 0.4478895
# [43] 0.4479795 0.4481595 0.4483395 0.4483395 0.4484295 0.4486095 0.4486095
# [50] 0.4486095 0.4486995 0.4490595 0.4493295 0.4495095 0.4497795 0.4497795
# [57] 0.4498695 0.4499595 0.4501395 0.4502295 0.4506795 0.4507695 0.4507695
# [64] 0.4507695 0.4510395 0.4513095 0.4515795 0.4520295 0.4521195 0.4521195
# [71] 0.4522095 0.4527495 0.4529295 0.4529295 0.4533795 0.4535595 0.4537395
# [78] 0.4542795 0.4543695 0.4549995 0.4551796 0.4553596 0.4557196 0.4558096
# [85] 0.4559896 0.4562596 0.4564396 0.4564396 0.4566196 0.4578796 0.4579696
# [92] 0.4581496 0.4582396 0.4583296 0.4585996 0.4589596 0.4594096 0.4596796
# [99] 0.4597696 0.4630996

 
#--------------------------------------------------
# task-2. random survival forest, get IBS error
#--------------------------------------------------

IBS.RSF.oob <- NULL
for(i in 1:100){
 rfsrc.obj <- rfsrc(surv.f, data=cr, nodesize=Nnode, ntree=Ntree, samptype="swr")
 ibs.time.all <- ibsRSF(rfsrc.obj)
 IBS.RSF.oob <- c(IBS.RSF.oob, ibs.time.all) 
}

# one typical run
# sort(IBS.RSF.oob)
#  [1] 0.1588891 0.1592320 0.1592588 0.1592786 0.1592820 0.1593686 0.1593687
#  [8] 0.1594684 0.1595015 0.1595018 0.1596188 0.1596371 0.1596594 0.1596712
# [15] 0.1596899 0.1596924 0.1597144 0.1597764 0.1597869 0.1598603 0.1598952
# [22] 0.1599071 0.1599207 0.1599367 0.1599626 0.1600347 0.1600406 0.1600499
# [29] 0.1600625 0.1600865 0.1600949 0.1601083 0.1601087 0.1601695 0.1601842
# [36] 0.1601878 0.1602107 0.1602223 0.1602249 0.1602452 0.1602475 0.1602685
# [43] 0.1603067 0.1603085 0.1603146 0.1603279 0.1603558 0.1603983 0.1604067
# [50] 0.1604168 0.1604477 0.1604485 0.1605178 0.1605502 0.1605841 0.1605959
# [57] 0.1606452 0.1606795 0.1606844 0.1606919 0.1607084 0.1607120 0.1607241
# [64] 0.1607602 0.1608451 0.1609154 0.1609162 0.1609472 0.1609811 0.1610024
# [71] 0.1610157 0.1610319 0.1610528 0.1610803 0.1611010 0.1611291 0.1611339
# [78] 0.1611764 0.1611958 0.1612090 0.1612167 0.1612175 0.1612273 0.1613046
# [85] 0.1613269 0.1613336 0.1614027 0.1614563 0.1615062 0.1615316 0.1615501
# [92] 0.1617066 0.1617376 0.1618476 0.1618646 0.1618810 0.1619498 0.1623804
# [99] 0.1624219 0.1625019


#---------------------------------------------------------------------
# task-3. Cox regression, get D-index error
# the "OOB" is created manually. remember we use ntree=500 in RSF, so here 
# we also repeat the Cox regression 500 times (which can be time-consuming)
#---------------------------------------------------------------------

D.COX.oob <- NULL
for(i in 1:100){

 di_tmp <- NULL
 # for(j in 1:500){
 for(j in 1:5){
  # note one sample can appear more than once in "train", because "with replacement"
  train <- sample(1:nrow(cr), nrow(cr), replace=TRUE)
  cox.obj <- coxph(surv.f, cr[train, ])

  # this is equivalent to cr[-unique(train),]
  # these will be "oob" samples in this j-iteration
  cr.val <- cr[-train, ] 
  newy <- Surv(cr.val$time, cr.val$status)
  newx <- subset(cr.val, select=-c(time, status))

  di_tmp <- c(di_tmp, dindexCOX(cox.obj, newx, newy)$dindex )
 }

 D.COX.oob <- c(D.COX.oob, mean(di_tmp) )
}

# one typical run
# sort(D.COX.oob)
#  [1] 0.4530422 0.4535487 0.4538870 0.4543685 0.4551308 0.4553555 0.4554340
#  [8] 0.4561047 0.4561185 0.4561227 0.4561803 0.4562450 0.4562584 0.4563706
# [15] 0.4565045 0.4565142 0.4565422 0.4566814 0.4567062 0.4568239 0.4568378
# [22] 0.4568984 0.4569381 0.4569755 0.4569834 0.4570843 0.4571400 0.4571824
# [29] 0.4572498 0.4572768 0.4573356 0.4575497 0.4576002 0.4576210 0.4576508
# [36] 0.4576822 0.4576867 0.4577025 0.4577452 0.4577459 0.4577753 0.4577783
# [43] 0.4578301 0.4579394 0.4579965 0.4580674 0.4580771 0.4581310 0.4581608
# [50] 0.4581721 0.4582025 0.4583051 0.4583885 0.4584062 0.4584664 0.4585556
# [57] 0.4586113 0.4586438 0.4586652 0.4588330 0.4590245 0.4590868 0.4591178
# [64] 0.4591448 0.4592044 0.4593708 0.4594110 0.4594338 0.4594682 0.4594755
# [71] 0.4595462 0.4595668 0.4596010 0.4596058 0.4596167 0.4596457 0.4598063
# [78] 0.4598413 0.4598826 0.4599302 0.4600689 0.4600890 0.4601392 0.4602334
# [85] 0.4602469 0.4602713 0.4603038 0.4605403 0.4606199 0.4606543 0.4606688
# [92] 0.4607049 0.4607301 0.4608225 0.4614698 0.4617318 0.4617749 0.4621043
# [99] 0.4626932 0.4632625


#---------------------------------------------------------------------
# task-4. Cox regression, get IBS error
# the "OOB" is created manually. remember we use ntree=500 in RSF, so here 
# we also repeat the Cox regression 500 times (which can be time-consuming)
#---------------------------------------------------------------------

IBS.COX.oob <- NULL
for(i in 1:100){

 ibs_tmp <- NULL
 # for(j in 1:500){
 for(j in 1:5){
  # note one sample can appear more than once in "train", because "with replacement"
  train <- sample(1:nrow(cr), nrow(cr), replace=TRUE)
  cox.obj <- coxph(surv.f, cr[train, ])

  # this is equivalent to cr[-unique(train),]
  # these will be "oob" samples in this j-iteration
  cr.val <- cr[-train, ]
  newy <- Surv(cr.val$time, cr.val$status)
  newx <- subset(cr.val, select=-c(time, status))
  ibs_tmp <- c(ibs_tmp,  ibsCOX(cox.obj, newx, newy) )
 }
 IBS.COX.oob <- c(IBS.COX.oob, mean(ibs_tmp))
}

# one typical run
# sort(IBS.COX.oob)
#  [1] 0.1966020 0.1966251 0.1966621 0.1966815 0.1969334 0.1973673 0.1974351
#  [8] 0.1974965 0.1974982 0.1975175 0.1975416 0.1975598 0.1976045 0.1976872
# [15] 0.1977783 0.1977786 0.1978179 0.1978341 0.1978757 0.1979029 0.1979618
# [22] 0.1979878 0.1979881 0.1980382 0.1980464 0.1980916 0.1980921 0.1981093
# [29] 0.1981977 0.1982144 0.1983032 0.1983189 0.1983556 0.1983830 0.1983996
# [36] 0.1984184 0.1984235 0.1984277 0.1984982 0.1985102 0.1985508 0.1985581
# [43] 0.1986052 0.1986149 0.1986167 0.1986231 0.1986359 0.1986469 0.1986888
# [50] 0.1986914 0.1987082 0.1987186 0.1987251 0.1987411 0.1987419 0.1987861
# [57] 0.1988932 0.1989593 0.1989986 0.1990140 0.1990333 0.1990421 0.1990824
# [64] 0.1991223 0.1991336 0.1991557 0.1991619 0.1992350 0.1992452 0.1992651
# [71] 0.1993178 0.1993197 0.1993271 0.1993376 0.1993438 0.1993595 0.1993667
# [78] 0.1993675 0.1993756 0.1994604 0.1994618 0.1994934 0.1995463 0.1995725
# [85] 0.1995855 0.1997427 0.1997870 0.1998226 0.1998343 0.1999007 0.1999464
# [92] 0.2000126 0.2001051 0.2001333 0.2002583 0.2002611 0.2003198 0.2003516
# [99] 0.2004234 0.2006202


