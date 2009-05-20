\name{sampalgonthefly-class}
\alias{sampalgonthefly}
\docType{class}
\alias{sampalgonthefly-class}

\title{Class "sampalgonthefly"}
\description{A sequential sampling algorithm that creates its boundaries
  on the fly.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{getalgonthefly(level,epsilon,halfspend)}.
}
\section{Slots}{
	 \describe{
    \item{\code{internal}:}{Object of class
      \code{"environment"}. Internal state of the algorithm. Do not access. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{sampalg}"}, directly.
}
\section{Methods}{
  \describe{
    \item{run}{\code{signature(alg = "sampalgonthefly")}: ... }
	 }
  \describe{
    \item{getboundaryandprob}{\code{signature(alg = "sampalgonthefly")}: ... }
	 }
}
\references{
  Gandy, A. (2009) Sequential Implementation of Monte Carlo Tests 
  with Uniformly Bounded Resampling Risk. To appear in JASA.
}
\author{ Axel Gandy }

\seealso{
 \code{\linkS4class{sampalgPrecomp}} 
}
\examples{
showClass("sampalgonthefly")
}
\keyword{classes}
