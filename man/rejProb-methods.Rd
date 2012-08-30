\name{rejProb-methods}
\alias{rejProb}
\docType{methods}
\alias{rejProb-methods}
\alias{rejProb,mmctestres-method}
\alias{rejProb,mmctest-method}
\title{Methods for class `mmctestres' and `mmctest', Package `simctest'}
\description{
  Function which returns empirical rejection probabilities.
  Threshold against e.g. 0.5 to obtain rejections (all rejProb>0.5 are rejected).
  Important: For usage in connection with thompson=TRUE (see the mmctest constructor).
}
\usage{
 rejProb(obj)
}
\arguments{
  \item{obj}{object of type `mmctestres' or `mmctest'.
  }
}
\section{Methods}{
  \describe{
    \item{rejProb(obj)}{ works with object of type mmctestres or mmctest. }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n[i])<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
  rejProb(a);
}
\keyword{methods}
