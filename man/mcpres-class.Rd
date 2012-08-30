\name{mcpres-class}
\docType{class}
\alias{mcpres-class}
\alias{show,mcpres-method}

\title{Class "mcpres"}
\description{Result returned by \code{mcp}}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("mcpres", ...)}.
}
\section{Slots}{
  \describe{       
    \item{\code{int}:}{Object of class \code{"numeric"}}
    \item{\code{cp}:}{Object of class \code{"numeric"}}
    \item{\code{beta}:}{Object of class \code{"numeric"}}
    \item{\code{N}:}{Object of class \code{"numeric"}}
    \item{\code{effort}:}{Object of class \code{"numeric"}}
    \item{\code{rescount}:}{Object of class \code{"numeric"}}
    \item{\code{truncated}:}{Object of class \code{"logical"}}
    \item{\code{taccepted}:}{Object of class \code{"logical"}}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "mcpres")}: ... }
	 }
}
\references{
  Gandy, A. and Rubin-Delanchy, P. An algorithm to compute the power of
  Monte Carlo tests with guaranteed precision. 	arXiv:1110.1248v1}
\author{Axel Gandy and Patrick Rubin-Delanchy}

\examples{
showClass("mcpres")
}
\keyword{classes}
