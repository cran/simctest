\name{sampalgPrecomp-class}
\alias{sampalgPrecomp}
\docType{class}
\alias{sampalgPrecomp-class}

\title{Class "sampalgPrecomp"}
\description{A sampling algorithm that precomputes the boundaries}
\section{Objects from the Class}{
Objects can be created by calls to \code{\link{getalgprecomp}}
}
\section{Slots}{
	 \describe{
    \item{\code{internal}:}{internal state of the object. Do not access. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{sampalg}"}, directly.
}
\references{
  Gandy, A. (2009) Sequential Implementation of Monte Carlo Tests 
  with Uniformly Bounded Resampling Risk. To appear in JASA.
}
\author{ Axel Gandy }


\examples{
showClass("sampalgPrecomp")
}
\keyword{classes}
