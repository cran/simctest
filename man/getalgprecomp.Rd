\name{getalgprecomp}
\alias{getalgprecomp}
\alias{getalgonthefly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct algorithms  ~~function to do ... ~~ }
\description{
  Constructs classes of type \code{\linkS4class{sampalgonthefly}} and \code{\linkS4class{sampalgPrecomp}}.
}
\usage{
getalgonthefly(level = 0.05, epsilon = 0.001, halfspend = 1000)
getalgprecomp(level = 0.05, epsilon = 0.001, halfspend = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{level}{the threshold.}
  \item{epsilon}{the bound on the resampling risk.}
  \item{halfspend}{number of steps after which half the error has been spent.}
}
\value{
  \code{getalgonthefly} returns an object of type \code{\linkS4class{sampalgonthefly}}.
  \code{getalgprecomp} returns an object of type \code{\linkS4class{sampalgPrecomp}}.
}
\references{
  Gandy, A. (2009) Sequential Implementation of Monte Carlo Tests 
  with Uniformly Bounded Resampling Risk. To appear in JASA.
}
\author{Axel Gandy}
\examples{
alg<-getalgprecomp()
run(alg, function() runif(1)<0.01)

alg<-getalgonthefly()
run(alg, function() runif(1)<0.01)
}
\keyword{classes}
