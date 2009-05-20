\name{getU-methods}
\alias{getU}
\docType{methods}
\alias{getU-methods}
\alias{getU,sampalgPrecomp-method}
\title{Methods for Function getU in Package `simctest'}
\description{
Returns the upper boundary for the stopping rule
}
\usage{
getU(alg,ind)
}
\arguments{
  \item{alg}{the sampling algorithm}
  \item{ind}{a vector of indices at which the upper stopping boundary should
    be returned}
}
\section{Methods}{
\describe{
\item{alg = "sampalgPrecomp"}{ ~~describe this method here }
}}
\examples{
   getU(getalgprecomp(),1:100)
}
\keyword{methods}
