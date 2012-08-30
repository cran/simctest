\name{mkdeltamid}
\alias{mkdeltamid}
\title{Function mkdeltamid in Package `simctest'}
\description{Easy creation of adaptive delta function}
\usage{
mkdeltamid(mindelta=0.02, maxdelta=0.1, llim=0.05, rlim=0.95)
}
\arguments{
  \item{mindelta}{desired length of CI for regions of interest, such as
    when the power is less than 0.05 or greater than 0.95.}
  \item{maxdelta}{desired length of CI when power is not in reregion of
    interest, e.g. between 0.05 and 0.95 }
  \item{llim}{change if want different left limit (i.e. not 0.05)}
  \item{rlim}{change if want different right limit (i.e. not 0.95)}
}

\value{
  A function, say \code{deltamid}, that specifies
  the user's desired precision depending on the midpoint of the computed confidence interval.
  If the current confidence interval has a midpoint M, then the
  algorithm will stop if deltamid(M) <= length of CI.
}
\references{
  Gandy, A. and Rubin-Delanchy, P. An algorithm to compute the power of
  Monte Carlo tests with guaranteed precision. arXiv:1110.1248v1}
\author{ Axel Gandy and Patrick Rubin-Delanchy}

\examples{
## only care about powers around 0.9 or higher
## (e.g. if want to check that the test is powerful enough).

deltamid <- mkdeltamid(mindelta=0.02, maxdelta=1, llim=0, rlim=0.9)

genstream <- function(){p <- runif(1); function(N){runif(N) <= p}}

## The power is 0.05. The algorithm should stop as soon as it is clear
## that the power is not larger than 0.9. (Must specify epsilon
## if using non-standard delta.)

res <- mcp(genstream, alpha=0.05, delta="adaptive", cp=0.99,
options=list(deltamid = deltamid, epsilon = 0.0001))

##should stop early.
res
}

