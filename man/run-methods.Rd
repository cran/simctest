\name{run-methods}
\alias{run}
\docType{methods}
\alias{run-methods}
\alias{run,sampalgPrecomp,ANY-method}
\alias{run,sampalgonthefly,ANY-method}
\alias{run,mmctest,mmctSamplerGeneric-method}

\title{Methods for Function run in Package `simctest'}
\description{Starts a sampling algorithm
}
\usage{
  run(alg,gensample,maxsteps)
}
\arguments{
  \item{alg}{the sampling algorithm. An object of type "sampalg" or "mmctest".
  }
  \item{gensample}{a function returing the result of one resampling step
    (0=no rejection, 1=rejection of the null hypothesis), or an object of
    type "mmctSamplerGeneric" if alg="mmctest".
  }
  \item{maxsteps}{the maximal number of steps to take
  }
}
\section{Methods}{
\describe{

\item{alg = "sampalgPrecomp"}{ the algorithm to be used }

\item{alg = "sampalgonthefly"}{ the algorithm to be used }

\item{alg = "mmctest", gensample = "mmctSamplerGeneric"}{ the algorithm to be used }
}}
\examples{
alg<-getalgonthefly()
res <- run(alg, function() runif(1)<0.2);
res
}
\keyword{methods}
