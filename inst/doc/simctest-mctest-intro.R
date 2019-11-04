### R code from vignette source 'simctest-mctest-intro.Rnw'

###################################################
### code chunk number 1: simctest-mctest-intro.Rnw:9-10
###################################################
options(width=80)


###################################################
### code chunk number 2: simctest-mctest-intro.Rnw:40-41
###################################################
library(simctest)


###################################################
### code chunk number 3: simctest-mctest-intro.Rnw:45-46 (eval = FALSE)
###################################################
## vignette("simctest-mctest-intro")


###################################################
### code chunk number 4: simctest-mctest-intro.Rnw:93-94
###################################################
gen <- function() { runif(1)<0.04 }


###################################################
### code chunk number 5: simctest-mctest-intro.Rnw:102-108
###################################################
J <- matrix(nrow=2,
            c(0,   1e-3,
              1e-3,1e-2,
              1e-2,0.05,
              0.05,1))
colnames(J) <- c("***","**","*","")


###################################################
### code chunk number 6: simctest-mctest-intro.Rnw:115-117
###################################################
res <- mctest(gen,J=Jstar,epsilon=0.001,batch=10,
	      batchincrement=1.1,maxbatch=100,method="simctest")


###################################################
### code chunk number 7: simctest-mctest-intro.Rnw:137-139
###################################################
mctest.RL(gen,J=Jstar,epsilon=0.001,batch=10,
	  batchincrement=1.1,maxbatch=100)


###################################################
### code chunk number 8: simctest-mctest-intro.Rnw:142-144
###################################################
mctest.simctest(gen,J=Jstar,epsilon=0.001,batch=10,
		batchincrement=1.1,maxbatch=100)


###################################################
### code chunk number 9: simctest-mctest-intro.Rnw:155-156
###################################################
res


###################################################
### code chunk number 10: simctest-mctest-intro.Rnw:167-172
###################################################
res$decision.interval
res$decision
res$est.p
res$batchedSamples
res$actualSamples


###################################################
### code chunk number 11: simctest-mctest-intro.Rnw:184-196
###################################################
dat <- matrix(nrow=5,ncol=7,byrow=TRUE,
  c(1,2,2,1,1,0,1,
    2,0,0,2,3,0,0,
    0,1,1,1,2,7,3,
    1,1,2,0,0,0,1,
    0,1,1,1,1,0,0))
loglikrat <- function(data) {
  cs <- colSums(data)
  rs <- rowSums(data)
  mu <- outer(rs,cs)/sum(rs)
  2*sum(ifelse(data<=0.5, 0,data*log(data/mu)))
}


###################################################
### code chunk number 12: simctest-mctest-intro.Rnw:201-208
###################################################
resample <- function(data){
  cs <- colSums(data)
  rs <- rowSums(data)
  n <- sum(rs)
  mu <- outer(rs,cs)/n/n
  matrix(rmultinom(1,n,c(mu)),nrow=dim(data)[1],ncol=dim(data)[2])
}


###################################################
### code chunk number 13: simctest-mctest-intro.Rnw:212-213
###################################################
t <- loglikrat(dat)


###################################################
### code chunk number 14: simctest-mctest-intro.Rnw:219-220
###################################################
gen <- function(){loglikrat(resample(dat))>=t}


###################################################
### code chunk number 15: simctest-mctest-intro.Rnw:224-227
###################################################
res <- mctest(gen,method="simctest")
mctest.simctest(gen)
mctest.RL(gen)


###################################################
### code chunk number 16: simctest-mctest-intro.Rnw:233-235
###################################################
res$decision.interval
res$decision


