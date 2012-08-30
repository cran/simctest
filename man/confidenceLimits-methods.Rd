\name{confidenceLimits-methods}
\alias{confidenceLimits}
\docType{methods}
\alias{confidenceLimits-methods}
\alias{confidenceLimits,mmctestres-method}
\alias{confidenceLimits,mmctest-method}
\title{Methods for class `mmctestres' and `mmctest', Package `simctest'}
\description{
  Function which returns a list containing lower confidence limits (vector `lowerLimits') and
upper confidence limits (vector `upperLimits').
}
\usage{
 confidenceLimits(obj)
}
\arguments{
  \item{obj}{object of type `mmctestres' or `mmctest'.
  }
}
\section{Methods}{
  \describe{
    \item{confidenceLimits(obj)}{ works with object of type mmctestres or mmctest. }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n)<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
  res <- confidenceLimits(a);
  lower <- res$lowerLimits;
  upper <- res$upperLimits;
}
\keyword{methods}
