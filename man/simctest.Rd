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
\seealso{}

\examples{
  #Example from Davison&Hinkley section 11.4, p. 534
  data(fir,package="boot")
  fir.mle <- c(sum(fir$count), nrow(fir))
  fir.gen <- function(data, mle){
    d<-data
    y <- sample(x=mle[2],size=mle[1],replace=TRUE)
    d$count <- tabulate(y,mle[2])
    d
  }
  fir.fun <- function(data)
    (nrow(data)-1)*var(data$count)/mean(data$count)
  obs<-fir.fun(fir)
  resampl<-function(){obs<fir.fun(fir.gen(data=fir,mle=fir.mle))}
  simctest(resampl)

  #now trying simctest.cont
  res <- simctest(resampl,maxsteps=10)
  res

  res <- cont(res,10)
  cont(res,10)

}
\author{Axel Gandy}
\keyword{}
