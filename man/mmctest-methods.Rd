\name{mmctest-methods}
\alias{mmctest}
\docType{methods}
\alias{mmctest-methods}

\title{Methods for class `mmctest', Package `simctest'}
\description{
  Constructor for class `mmctest'.
}
\usage{
 mmctest(epsilon=0.01, threshold=0.1, h)
}
\arguments{
  \item{epsilon}{threshold for testing.
  }
  \item{threshold}{threshold for testing.
  }
  \item{h}{reference to a multiple testing function of the
form function(p, threshold) which returns the set of
rejected indices.
  }
}
\section{Methods}{
  \describe{
    \item{mmctest(epsilon=0.01, threshold=0.1, h)}{ returns object of type `mmctest'. }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n)<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
}
\keyword{methods}
