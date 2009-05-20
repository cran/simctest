\name{sampalgontheflyres-class}
\docType{class}
\alias{sampalgontheflyres-class}
\alias{contalg,sampalgontheflyres-method}

\title{Class "sampalgontheflyres"}
\description{Class returned as result from \code{simctest} and \code{run}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sampalgontheflyres", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{porig}:}{Object of class \code{"numeric"}  }
    \item{\code{U}:}{Object of class \code{"numeric"}  }
    \item{\code{L}:}{Object of class \code{"numeric"}  }
    \item{\code{ind}:}{Object of class \code{"numeric"}  }
    \item{\code{preverr}:}{Object of class \code{"numeric"}  }
    \item{\code{p.value}:}{Object of class \code{"numeric"}  }
    \item{\code{steps}:}{Object of class \code{"numeric"}  }
    \item{\code{pos}:}{Object of class \code{"numeric"}  }
    \item{\code{alg}:}{Object of class \code{"sampalg"}  }
    \item{\code{gen}:}{Object of class \code{"function"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{sampalgres}"}, directly.
}
\section{Methods}{
  \describe{
    \item{contalg}{\code{signature(data = "sampalgontheflyres")}: ... }
	 }
}
\references{
  Gandy, A. (2009) Sequential Implementation of Monte Carlo Tests 
  with Uniformly Bounded Resampling Risk. To appear in JASA.
}
\author{ Axel Gandy }
\note{ ~~further notes~~ }


\seealso{
  \code{\link{simctest}}, \code{\linkS4class{sampalgres}} 
}
\examples{
showClass("sampalgontheflyres")
}
\keyword{classes}
