\name{mmctSampler-class}
\docType{class}
\alias{mmctSampler-class}

\title{Class "mmctest"}
\description{Wrapper-Class for "mmctestInterfaceGeneric", takes a function, the number
of hypotheses and returns derived object of class "mmctestInterfaceGeneric".
Class provides a slot for additional data.
The function f(ind,n,data) has to return n[i] new samples for each hypothesis ind[i]
in vector "ind", where i=1...length(ind).
The data stored in the data slot of class "mmctSampler" is also
passed on to "f".}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{mmctSampler(f=...,num=...,data=...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{f}:}{Object of class \code{"function"}  }
    \item{\code{num}:}{Object of class \code{"numeric"}  }
    \item{\code{data}:}{Object of class \code{"numeric"}  }
  }
}
\section{Methods}{
  \describe{
    \item{getSamples}{\code{signature(obj="mmctSampler", ind="numeric", n="numeric")}: ... }
    \item{getNumber}{\code{signature(obj="mmctSampler")}: ... }
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
}
\keyword{classes}
