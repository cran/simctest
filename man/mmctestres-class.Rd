\name{mmctestres-class}
\docType{class}
\alias{mmctestres-class}
\alias{show,mmctestres-method}

\title{Class "mmctestres"}
\description{Class which stores current result of type "mmctest".}
\section{Objects from the Class}{
Objects should not be created directly.
Objects returned by calls of the form \code{new("mmctest", ...)}
are of type mmctestres.
}
\section{Slots}{
	 \describe{
    \item{\code{internal}:}{Object of class \code{"environment"}  }
    \item{\code{epsilon}:}{Object of class \code{"numeric"}  }
    \item{\code{threshold}:}{Object of class \code{"numeric"}  }
    \item{\code{r}:}{Object of class \code{"numeric"}  }
    \item{\code{R}:}{Object of class \code{"numeric"}  }
    \item{\code{h}:}{Object of class \code{"function"}  }
    \item{\code{gensample}:}{Object of class \code{"mmctSamplerGeneric"}  }
    \item{\code{g}:}{Object of class \code{"numeric"}  }
    \item{\code{num}:}{Object of class \code{"numeric"}  }
    \item{\code{A}:}{Object of class \code{"numeric"}  }
    \item{\code{B}:}{Object of class \code{"numeric"}  }
    \item{\code{C}:}{Object of class \code{"numeric"}  }
    \item{\code{thompson}:}{Object of class \code{"logical"}  }
    \item{\code{rejprob}:}{Object of class \code{"logical"}  }
  }
}
\section{Methods}{
  \describe{
    \item{mainalg}{\code{signature(obj = "mmctestres", stopcrit = "numeric")}: ... }
    \item{cont}{\code{signature(data = "mmctestres", steps = "numeric")}: ... }
    \item{show}{\code{signature(object = "mmctestres")}: ... }
    \item{pEstimate}{\code{signature(obj = "mmctestres")}: ... }
    \item{rejProb}{\code{signature(obj = "mmctestres")}: ... }
    \item{confidenceLimits}{\code{signature(obj = "mmctestres")}: ... }
    \item{testResult}{\code{signature(obj = "mmctestres")}: ... }
    \item{summary.mmctestres}{\code{signature(object = "mmctestres")}: ... }
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
  # a is object of type "mmctestres" now
}
\keyword{classes}
