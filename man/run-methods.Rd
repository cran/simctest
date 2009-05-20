\name{run-methods}
\alias{run}
\docType{methods}
\alias{run-methods}
\alias{run,sampalgPrecomp-method}
\alias{run,sampalgonthefly-method}
\title{Methods for Function run in Package `simctest'}
\description{Starts a sampling algorithm
}
\usage{
  run(alg,gensample,maxsteps=1e4)
}
\arguments{
  \item{alg}{the sampling algorithm. An object of type sampalg.
  }
  \item{gensample}{a function returing the result of one resampling step
    (0=no rejection, 1=rejection of the null hypothesis)
  }
  \item{maxsteps}{the maximal number of steps to take
  }
}
\section{Methods}{
\describe{

\item{alg = "sampalgPrecomp"}{ ~~describe this method here }

\item{alg = "sampalgonthefly"}{ ~~describe this method here }
}}
\examples{
alg<-getalgonthefly()
res <- run(alg, function() runif(1)<0.2);
res
}
\keyword{methods}
