\name{sampalgres-class}
\docType{class}
\alias{sampalgres-class}
\alias{contalg,sampalgres-method}
\alias{show,sampalgres-method}

\title{Class "sampalgres"}
\description{Results returned by \code{run} - Internal.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sampalgres", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{p.value}:}{Object of class \code{"numeric"}  }
    \item{\code{steps}:}{Object of class \code{"numeric"}  }
    \item{\code{pos}:}{Object of class \code{"numeric"}  }
    \item{\code{alg}:}{Object of class \code{"sampalg"}  }
    \item{\code{gen}:}{Object of class \code{"function"}  }
  }
}
\section{Methods}{
  \describe{
    \item{confint}{\code{signature(object = "sampalgres", parm = "missing")}: ... }
    \item{contalg}{\code{signature(data = "sampalgres")}: ... }
    \item{getbounds}{\code{signature(data = "sampalgres")}: ... }
    \item{show}{\code{signature(object = "sampalgres")}: ... }
	 }
}
\references{
  Gandy, A. (2009) Sequential Implementation of Monte Carlo Tests 
  with Uniformly Bounded Resampling Risk. To appear in JASA.}
\author{ Axel Gandy }

\examples{
showClass("sampalgres")
}
\keyword{classes}
