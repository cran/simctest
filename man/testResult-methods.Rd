\name{testResult-methods}
\alias{testResult}
\docType{methods}
\alias{testResult-methods}
\alias{testResult,mmctestres-method}
\alias{testResult,mmctest-method}
\title{Methods for class `mmctestres' and `mmctest', Package `simctest'}
\description{
  Function which returns a list containing indices of rejected hypotheses (vector `rejected'),
nonrejected hypotheses (vector `nonrejected') and undecided hypotheses (vector `undecided').
}
\usage{
 testResult(obj)
}
\arguments{
  \item{obj}{object of type `mmctestres' or `mmctest'.
  }
}
\section{Methods}{
  \describe{
    \item{testResult(obj)}{ works with object of type mmctestres or mmctest. }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n)<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
  res <- testResult(a);
  rejected <- res$rejected;
  nonrejected <- res$nonrejected;
  undecided <- res$undecided;
}
\keyword{methods}
