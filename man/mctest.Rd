\name{mctest}
\alias{mctest}
\alias{mctest.RL}
\alias{mctest.simctest}
\alias{print.mctestres}
\alias{J}
\alias{Jstar}
\alias{Jnarrow}
\title{Sequential implementation of Monte Carlo
  tests with multiple thresholds}
\description{
  Sequential implementation
  of the Monte Carlo test with multiple thresholds.

  Implementation of the Robbins-Lai (mctest.RL) and SIMCTEST (mctest.simctest) approaches to compute a decision interval (and decision) with respect to several thresholds/ p-value buckets.
  The function "mctest" is a wrapper function for both the Robbins-Lai and the SIMCTEST approach which calls one of the two using an additional parameter "method" (method="simctest" for SIMCTEST and method="RL" for Robbins-Lai).
}
\usage{
mctest(gen,J=Jnarrow,epsilon=0.001,batch=10,batchincrement=1.1,maxbatch=100,
  method=c("simctest","RL"))
mctest.RL(gen,J=Jnarrow,epsilon=0.001,batch=10,batchincrement=1.1,maxbatch=100)
mctest.simctest(gen,J=Jnarrow,epsilon=0.001,batch=10,batchincrement=1.1,maxbatch=100)
J
Jstar
Jnarrow
\method{print}{mctestres}(x,...)
}
\arguments{
  \item{gen}{function that performs one sampling step. 
    Returns 0 (sampled test statistic does not exceed the observation)
    or 1 (sampled test static exceeds the observation)
  }
  \item{method}{which method to use for stopping}
  \item{J}{interval buckets to use. A matrix with two rows, each column
    describes an interval bucket. Column names give the code for the
    interval bucket. Defaults to Jnarrow
  }
  \item{epsilon}{error bound}
  \item{batch}{initial number of samples to use before checking for stopping}
  \item{batchincrement}{factor by which the batch size gets multiplied
    after each step. 1 would mean no increment}
  \item{maxbatch}{maximum batch size}
  \item{x}{object of type "mctestres"}
  \item{...}{further arguments}
}
\value{
     \code{mctest}, \code{mctest.RL} and \code{mctest.simctest} all
  return an object of class type  \code{mctestres}, which has a print
  function (\code{print.mctestres}).

   An object of class \code{mctestres} is a list with  the
   following components: \code{step} (total batched number of samples drawn), \code{decision.interval} (interval for the p-value), \code{decision} (expressing significance), \code{est.p} (an estimate of the p-value) and \code{realn} (the actual number of samples taken without batching).
}
\references{
  Ding, D., Gandy, A. and Hahn, G. (2017) Implementing Monte Carlo Tests
  with Multiple Thresholds. Working paper.
}
\examples{
  #Example used in the above paper
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
  mctest(gen)
  mctest.simctest(gen)
  mctest.RL(gen)
}
