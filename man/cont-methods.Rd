\name{cont-methods}
\alias{cont}
\docType{methods}
\alias{cont-methods}
\alias{cont,sampalgres-method}
\alias{cont,sampalgontheflyres-method}
\title{Methods for Function `cont' in Package `simctest'}
\description{
  Continues the sampling for some more steps.
}
\usage{
 cont(data,steps)
}
\arguments{
  \item{data}{a result of a run of a sampling algorithm that has not
    come to a conclusion yet.
  }
  \item{steps}{maximum number of further iterations to take.
  }
}
\section{Methods}{
  \describe{
    \item{data = "sampalgres"}{ works with the algorithm based on precomputation. }   
    \item{data = "sampalgontheflyres"}{ works with the on-the-fly algorithm. }
  }
}
\examples{
res <- simctest(function() runif(1)>0.95,maxsteps=10);
res
res <- cont(res,1000)
res
res <- cont(res,1000)
res
}
\keyword{methods}
