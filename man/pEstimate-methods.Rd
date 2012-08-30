\name{pEstimate-methods}
\alias{pEstimate}
\docType{methods}
\alias{pEstimate-methods}
\alias{pEstimate,mmctestres-method}
\alias{pEstimate,mmctest-method}
\title{Methods for class `mmctestres' and `mmctest', Package `simctest'}
\description{
  Function which shows current estimates of p-values.
}
\usage{
 pEstimate(obj)
}
\arguments{
  \item{obj}{object of type `mmctestres' or `mmctest'.
  }
}
\section{Methods}{
  \describe{
    \item{pEstimate(obj)}{ works with object of type mmctestres or mmctest. }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n[i])<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
  pEstimate(a);
}
\keyword{methods}
