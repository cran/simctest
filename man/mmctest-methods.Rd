\name{mmctest-methods}
\alias{mmctest}
\docType{methods}
\alias{mmctest-methods}

\title{Methods for class `mmctest', Package `simctest'}
\description{
  Constructor for class `mmctest'.
}
\usage{
 mmctest(epsilon=0.01, threshold=0.1, r=10000, h, thompson=F, R=1000)
}
\arguments{
  \item{epsilon}{probability of any misclassification one is willing to tolerate
  }
  \item{threshold}{threshold for testing.
  }
  \item{r}{parameter of the spending sequence, see vignette
  }
  \item{h}{reference to a multiple testing function of the
form function(p, threshold) which returns the set of
rejected indices.
  }
  \item{thompson}{if set to true, mmctest will use a Thompson strategy to draw further samples
  }
  \item{R}{number of repetitions (=draws from the posterior distributions) used to calculate
  empirical probabilities of each hypothesis being rejected -- used to calculate weights
  in QuickMMCTest (option thompson=TRUE in the mmctest constructor)
  }
}
\section{Methods}{
  \describe{
    \item{mmctest(epsilon=0.01, threshold=0.1, r=10000, h)}{ returns object of type `mmctest'. }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n[i])<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
}
\keyword{methods}
