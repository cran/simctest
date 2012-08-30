### R code from vignette source 'simctest-intro.Rnw'

###################################################
### code chunk number 1: simctest-intro.Rnw:11-12
###################################################
options(width=80)


###################################################
### code chunk number 2: simctest-intro.Rnw:69-70
###################################################
library(simctest)


###################################################
### code chunk number 3: simctest-intro.Rnw:74-75 (eval = FALSE)
###################################################
## vignette("simctest-intro")


###################################################
### code chunk number 4: simctest-intro.Rnw:92-94
###################################################
res <- simctest(function() runif(1)<0.07);
res


###################################################
### code chunk number 5: simctest-intro.Rnw:98-99
###################################################
confint(res)


###################################################
### code chunk number 6: simctest-intro.Rnw:107-109
###################################################
res <- simctest(function() runif(1)<0.05);
res


###################################################
### code chunk number 7: simctest-intro.Rnw:115-117
###################################################
res <- cont(res,10000)
res


###################################################
### code chunk number 8: simctest-intro.Rnw:124-152
###################################################
  dat <- matrix(nrow=5,ncol=7,byrow=TRUE,
                c(1,2,2,1,1,0,1, 2,0,0,2,3,0,0, 0,1,1,1,2,7,3, 1,1,2,0,0,0,1, 0,1,1,1,1,0,0))
  loglikrat <- function(data){
    cs <- colSums(data)
    rs <- rowSums(data)
    mu <- outer(rs,cs)/sum(rs)
    2*sum(ifelse(data<=0.5, 0,data*log(data/mu)))
  }
  resample <- function(data){
    cs <- colSums(data)
    rs <- rowSums(data)
    n <- sum(rs)
    mu <- outer(rs,cs)/n/n
    matrix(rmultinom(1,n,c(mu)),nrow=dim(data)[1],ncol=dim(data)[2])
  }
  t <- loglikrat(dat);

  # function to generate samples
  gen <- function(){loglikrat(resample(dat))>=t}

  #using simctest
  simctest(gen,maxsteps=10000)

  #now trying simctest.cont
  res <- simctest(gen,maxsteps=500)
  res

  cont(res,20000)


###################################################
### code chunk number 9: simctest-intro.Rnw:180-189
###################################################
genstream <- function(){
  D <- list(group1 = rnorm(8, mean=0), group2 = rnorm(4, mean=0.5))
  t <- mean(D$group2) - mean(D$group1)
  stream <- function(){
    group <- (c(D$group1, D$group2))[sample.int(12)]
    T <- mean(group[9:12]) - mean(group[1:8])
    T >= t      
  }    
}


###################################################
### code chunk number 10: simctest-intro.Rnw:206-207
###################################################
genstream <- function(){p <- runif(1); function(N){runif(N) <= p}}


###################################################
### code chunk number 11: simctest-intro.Rnw:214-216
###################################################
res<-mcp(genstream, options=list(reports=FALSE))
res


###################################################
### code chunk number 12: simctest-intro.Rnw:223-225
###################################################
res<-mcp(genstream, delta=0.05, options=list(reports=FALSE))
res


###################################################
### code chunk number 13: simctest-intro.Rnw:230-234
###################################################
res <- mcp(genstream, options=list(reports=FALSE, 
                        deltamid=mkdeltamid(mindelta=0.02, maxdelta=1, llim=0, rlim=0.9),
                        epsilon=0.0001))
res


