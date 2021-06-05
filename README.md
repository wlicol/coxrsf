
These are the R scripts used in this paper:

On fair performance comparison between random survival forest and Cox regression: an example of colorectal cancer study

Sirin Cetin, Ayse Ulgen, Isa Dede, Wentian Li

SciMedicine Journal, 3(1):66-76 (2021).

Data file: colorectal.txt
-------------------------

These information are provided for 221 colorectal cancer patients:

```
 1. age
 2. gender (1=male, 2=female)
 3. cancer type (1=colon, 2=rectum)
 4. leukocyte count
 5. neutrophile count
 6. lymphocyte count
 7. platelet count
 8. event(death) status (1=dead, 0=survived)
 9. time-to-event (days)
```

R script: IBS-and-D-for-RSF-and-Cox.R
-------------------------------------

There are two commonly used measures for survival analysis errors: discordance
index (D-index) which is 1 minus concordance index (1-Cindex), and integrated
Brier's score (IBS). 

If we use the "rfsrc" function from "randomForestSRC" R package, the resulting object would
contain the Dindex for out-of-bag (OOB) samples: rfsrc.obj$err.rate[rfsrc.obj$ntree] ,
and there is no need to provide extra script.

Our ibsRSF function calculates IBS error from a rfsrc object. There is an extra option (if
ave.flag=0) to show Brier's score at each time point. Another extra option (if self=1) is
to calculate the IBS for in-the-bag (training) samples.

If we use "coxph" function from the "survival" R package, our dindexCOX function calculates
Dindex for new samples (something like out-of-bag samples). An extra option (if self=1)
is to calculate Dindex for in-the-bad (training) samples. 

If we use "coxph" function from the "survival" R package, our ibsCOX function calculates
IBS for new samples (something like out-of-bag samples). There is an extra option (if
ave.flag=0) to show Brier's score at each time point. 

R function: ibsRSF
------------------

ibsRSF \<- function(rfsrc.obj, self=0, ave.flag=1)

arguments:

```
 rfsrc.obj: object from a rfsrc() run

 ave.flag: 1 (default) if the output is the average of Brier's score at all time points, 0 if the output is the time
series of Brier's score.

 self: 0 (default) if we use the out-of-bag samples to calculate IBS, 1 if we use the in-the-bag (training) samples.
```

R function: dindexCOX
---------------------

dindexCOX \<- function(cox.obj, newx, newy, self=0)

arguments:

```
 cox.obj: object from a coxph() run

 newx:  a matrix representing the independent variables (number of rows=number of samples, number of columns= number of
independent variables)

 newy: a two-column matrix representing the dependent variable (number of rows=number of samples, first column
is time-to-event data, second column is the event status)

 self: 0 (default) if we use the out-of-bag samples to calculate IBS, 1 if we use the in-the-bag (training) samples.



R function: ibsCOX
------------------

ibsCOX \<- function(cox.obj, newx, newy, ave.flag=1)

arguments:

```
  cox.obj: object from a coxph() run

 newx:  a matrix representing the independent variables (number of rows=number of samples, number of columns= number of
independent variables)

 newy: a two-column matrix representing the dependent variable (number of rows=number of samples, first column
is time-to-event data, second column is the event status)

 ave.flag: 1 (default) if the output is the average of Brier's score at all time points, 0 if the output is the time
series of Brier's score.
```


R script: example.R
-------------------

This R script illusrtates the use of the above three R functions by applying them to the colorectal survival data.


R script: some-calculation-used-in-the-paper.R
----------------------------------------------

This R script illustrates some calculation used in the Figures of our paper, include:

```
 tuning parameters (ntree, nodesize) in RSF

 "0.632" weighting of out-of-bag and in-the-bag IBS errors in RSF

 "0.632" weighting of out-of-bag and in-the-bag IBS errors in Cox regression

 10-fold IBS errors in Cox regression
```


