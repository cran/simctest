\name{simctest}
\alias{simctest}
\title{Sequential implementation of Monte Carlo
  tests}
\description{
  Wrapper function for convenient use of the sequential implementation
  of the Monte Carlo test.
}
\usage{
simctest(gensample, level=0.05, epsilon=1e-3, maxsteps=1e4)
}
\arguments{

  \item{gensample}{function that performs one sampling step. 
    Returns 0 (sampled test statistic does not exceed the observation)
    or 1 (sampled test static exceeds the observation).
  }
  \item{level}{level passed to \code{\link{getalgonthefly}}}
  \item{epsilon}{error bound epsilon passed to \code{\link{getalgonthefly}}}
  \item{maxsteps}{maximal number of steps to take}

}
\value{
  An object of class \code{\linkS4class{sampalgres}}.
}
\references{
  Gandy, A. (2009) Sequential Implementation of Monte Carlo Tests 
  with Uniformly Bounded Resampling Risk. To appear in JASA.
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
  simctest(gen,maxsteps=10000)

  #now trying simctest.cont
  res <- simctest(gen,maxsteps=500)
  res

  cont(res,20000)

}
\author{Axel Gandy}
