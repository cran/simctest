# multithresh.R
# implementation of "Implementing Monte Carlo Tests with Multiple Thresholds" with both Robbins-Lai and simctest confidence intervals
# file contains 2 classes

# *************************************
# TO EXPORT: J and Jstar
#            mctest (main routine)
#            mctest.RL, mctest.simctest
# *************************************

# definition of the thresholds
J <- matrix(nrow=2,
            c(0,   1e-3,
              1e-3,1e-2,
              1e-2,0.05,
              0.05,1))
colnames(J) <- c("***","**","*","")

Jstar <- matrix(nrow=2,
                c(0,   1e-3,
                  5e-4,2e-3,
                  1e-3,0.01,
                  8e-3,0.012,
                  0.01,0.05,
                  0.045,0.055,
                  0.05,1))
colnames(Jstar) <- c("***","**~","**","*~","*","~","")



# mctest

# multiple thresholds based on ROBBINS-LAI
mctest.RL <- function(gen,J=Jstar,epsilon=0.001,batch=10,batchincrement=1.1,maxbatch=100) {
  Sn = 0
  n = 0
  while(TRUE){
    n = n+batch
    Sn <- Sn + sum(replicate(batch,gen()))
    stopcrit <- (n+1)*dbinom(Sn,prob=J,size=n)<=epsilon
    deriv <- Sn/J-(n-Sn)/(1-J)
    if (any(J[1,]==0.)){
        stopcrit[1,J[1,]==0.] <- TRUE
        deriv[1,J[1,]==0.] <- 1
    }
    if (any(J[2,]==1.)){
        stopcrit[2,J[2,]==1.] <- TRUE
        deriv[2,J[2,]==1.] <- -1
    }
    decided <- stopcrit[1,]&stopcrit[2,]&deriv[1,]>=0&deriv[2,]<=0
    if (any(decided)){
        j <- min(which(decided))
        res <- list(batchedSamples=n,decision.interval=J[,j],decision=colnames(J)[j],est.p=Sn/n,actualSamples=n)
        class(res) <- "mctestres"
        return(res)
    }
    batch <- floor(min(maxbatch,batch*batchincrement))
  }
}
 
# multiple thresholds based on SIMCTEST
mctest.simctest <- function(gen,J=Jstar,epsilon=0.001,batch=10,batchincrement=1.1,maxbatch=100) {
  unique.alphas <- sort(unique(as.vector(J)))
  unique.alphas <- unique.alphas[-c(1,length(unique.alphas))]
  firststep <- TRUE
  I <- c(0,1)
  genloc <- function() {i <<- i+1;sample[i]}
  n <- 0 # number of steps before stopping
  realn <- 0
  while(TRUE){
      # run or continue all algorithms
      sample <- replicate(batch,gen())
      realn <- realn+batch
      if (firststep){
          allres <- list()
          for (j in 1:length(unique.alphas)){
              i <- 0;
              allres[[j]] <- simctest(genloc,epsilon=epsilon/2,level=unique.alphas[j],maxsteps=batch)
          }
      }else{
          for (j in 1:length(allres)){
              i <- 0;             
              allres[[j]] <- cont(allres[[j]],batch)
          }
      }
      #evaluate outcomes

      j <- 1
      while (j<=length(allres)){
          object <- allres[[j]]
          if (!is.na(object@p.value)){
              #stopped - adjust interval estimate
              if (object@p.value<=unique.alphas[j]){
                  #below threshold
                  I[2]<-min(I[2],unique.alphas[j])
              }else{
                  I[1]<-max(I[1],unique.alphas[j])
              }
              # remember number of steps taken
              if (allres[[j]]@steps>=n){
                  n <- allres[[j]]@steps
                  Sn <- allres[[j]]@pos
              }
              # now delete alpha from list
              allres[[j]] <- NULL   
              unique.alphas <- unique.alphas[-j]
          }else{
              j <- j+1
          }          
      }

      # now check if we are done
      decided <- J[1,]<=I[1] & I[2]<=J[2,]
      if (any(decided)){
          j <- min(which(decided))
          res <- list(batchedSamples=n,decision.interval=J[,j],decision=colnames(J)[j],est.p=Sn/n,actualSamples=realn)
          class(res) <- "mctestres"
          return(res)
      }
      #adjust batch length
      batch <- floor(min(maxbatch,batch*batchincrement))
      firststep <- FALSE
  }
}

# combined function
# returns computation as mctestres object
mctest <- function(gen,J=Jstar,epsilon=0.001,batch=10,batchincrement=1.1,maxbatch=100,method=c("simctest","RL")) {
  method <- match.arg(method)
  if(method=="simctest")  return(mctest.simctest(gen=gen,J=J,epsilon=epsilon,batch=batch,batchincrement=batchincrement,maxbatch=maxbatch))
  if(method=="RL")  return(mctest.RL(gen=gen,J=J,epsilon=epsilon,batch=batch,batchincrement=batchincrement,maxbatch=maxbatch))
  stop("Unknown method.")
}

print.mctestres <- function(x,...) {
  cat(paste("Interval for p-value: [",x$decision.interval[1],",",x$decision.interval[2],"]\n",sep=""))
  cat(paste("Decision: ",x$decision,"\n",sep=""))
  cat(paste("Estimate of p-value: ",x$est.p,"\n",sep=""))
  cat(paste("Batched number of samples: ",x$batchedSamples,"\n",sep=""))
  cat(paste("Actual number of samples: ",x$actualSamples,"\n",sep=""))
}
