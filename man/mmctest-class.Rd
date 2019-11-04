\name{mmctest-class}
\docType{class}
\alias{mmctest-class}

\title{Class "mmctest"}
\description{Class which creates an object of type "mmctestres".}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{mmctest(h=...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{internal}:}{Object of class \code{"environment"}  }
  }
}
\section{Methods}{
  \describe{
    \item{run}{\code{signature(alg = "mmctest", gensample = "mmctSamplerGeneric", maxsteps = "numeric")}: ... }
	 }
}
\references{
  Gandy, A. and Hahn, G. (2014) MMCTest - a safe algorithm for
  implementing multiple Monte Carlo tests. Scandinavian Journal of
  Statistics, 41(4):1083--1101
}
\author{ Axel Gandy and Georg Hahn }

\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n[i])<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
}
\keyword{classes}
