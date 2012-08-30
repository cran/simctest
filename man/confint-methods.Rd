\name{confint-methods}
\alias{confint}
\docType{methods}
\alias{confint-methods}
\alias{confint,ANY,ANY-method}
\alias{confint,sampalgres,missing-method}
\title{Methods for Function run in Package `simctest'}
\description{
Computes a confidence interval for the p-value
}
\usage{
confint(object,parm,level=0.95,...)
}
\arguments{
  \item{object}{An object of type \code{\linkS4class{sampalgres}}
    resulting from a previous call to \code{\link[=run-methods]{run}}
    or
    \code{\link[=cont-methods]{cont}}.
  }
  \item{parm}{must be missing.}
  \item{level}{the desired coverage probability.}
  \item{...}{additional argument(s). Currently not used}
}
\section{Methods}{
  \describe{
    
    \item{object = "ANY", parm = "ANY"}{  Generic function: see \code{\link{confint}}. }
    
    \item{object = "sampalgres", parm = "missing"}{Computes a confidence
      interval for the p-value with the coverage probability given by \code{level}. }
}}
\examples{
alg<-getalgonthefly()
res <- run(alg, function() runif(1)<0.05);
res
confint(res)
}
\keyword{methods}
