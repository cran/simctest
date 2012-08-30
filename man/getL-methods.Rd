\name{getL-methods}
\alias{getL}
\docType{methods}
\alias{getL-methods}
\alias{getL,sampalgPrecomp-method}
\title{Methods for Function getL in Package `simctest'}
\description{
Returns the lower boundary for the stopping rule
}
\usage{
##S4 method
getL(alg,ind)
}
\arguments{
  \item{alg}{the sampling algorithm}
  \item{ind}{a vector of indices at which the lower stopping boundary should
    be returned}
}
\section{Methods}{
\describe{
\item{alg = "sampalgPrecomp"}{ the sampling algorithm to be used }
}}
\examples{
   getL(getalgprecomp(),1:100)
}
\keyword{methods}
