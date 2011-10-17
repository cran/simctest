\name{mcp}
\alias{mcp}
\title{Function mcp in package `simctest'}
\description{An algorithm for the computation of the power of Monte Carlo tests with guaranteed precision}
\usage{
mcp(genstream,alpha=0.05,delta="adaptive",
cp=0.99,maxeffort=Inf,options = list())
}
\arguments{
  \item{genstream}{a function that returns a function that returns a
    random Bernoulli variable (each stream corresponds to a
    dataset. 0 = (T<t), 1= (T>=t) where t is computed from the dataset
    and T is a resampled test-statistic from that dataset.)}
  \item{alpha}{the level of the test.}
  \item{delta}{the desired length of confidence interval, or "adaptive" if using adaptive delta. See details.}
  \item{maxeffort}{maximum effort. Effort is total number of samples taken. Set to finite value if needed (the resulting confidence interval still has the guaranteed coverage probability, but may not be as `short' as desired). Can also interrupt the algorithm during main loop and get a result of class \code{"mcpres"}.}
  \item{cp}{the desired coverage probability.}
  \item{options}{Additional options. See details}
}
\details{
  \code{options$maxeffort}: set to maximum allowable effort.
  
  \code{options$reports}: set to \code{FALSE} if onscreen reports are not wanted. 
  
  \code{options$file}: optional file-name to save results to.
  
  \code{options$pilotn}: number of streams in pilot (1000 by default).
  
  \code{options$pilotmaxsteps}: maxsteps in pilot (1000 by default).
  
  \code{options$gammapilotprop}: proportion of error spent on pilot CI
  (0.1 by default)
  
  \code{options$gammatestprop}: proportion of error spent on testing
  remaining paths (default is 0.1)

  \code{options$spendgammatest}: spending sequence for the testing
    procedure on the remaining streams. Must be a non-negative function
    of integers with positive limit 1 (\eqn{t/(20+t)} by default).
    
  \code{options$eta}: internal parameter to the testing procedure on the
  remaining streams (0.05 by default).
  
  \code{options$maxstepsbase}: initial maximum number of steps (500 by
  default)

  \code{options$maxstepsinc}: multiplier for the maximum number of steps
  thereafter (1.5 by default).
  
  \code{options$maxbatch}: multiplier for the maximum number of steps
  thereafter (200000 by default).

  \code{options$deltamid}: adaptive delta function. Describes the length
  of the confidence interval desired depending on the midpoint of the
  interval. By default the function requires 0.02 for intervals
  containing 0.05 or lower or 0.95 or higher, and 0.1 otherwise. If
  using non-default adaptive delta must also specify epsilon
  (below).
  
  \code{options$epsilon}: error probability for each stream. Only set if using non-standard adaptive delta.
}

\seealso{\code{mkdeltamid}}

\value{
  An object of class \code{"mcpres"} with slots:
  \item{int}{confidence interval for power.}
  \item{cp}{coverage probability.}
  \item{beta}{Estimate of power.}
  \item{N}{the number of streams started in main loop (or in pilot if
    stopped after pilot).}
  \item{effort}{total number of samples generated.}
  \item{rescount}{number of positive and negative outcomes.}
  \item{truncated}{boolean indicating whether procedure was truncated by
  user-specified maxeffort.}
  \item{taccepted}{boolean indicating whether the procedure stopped as a
    result of a hypothesis test or brute force (the confidence interval coverage probability is guaranteed in either case.)}
}
\references{
  Gandy, A. and Rubin-Delanchy, P. An algorithm to compute the power of
  Monte Carlo tests with guaranteed precision. arXiv:1110.1248v1}
\author{Axel Gandy and Patrick Rubin-Delanchy}

\examples{
#Example where we know the power should be the level of the test
genstream <- function(){p <- runif(1); function(N){runif(N) <= p}}

res <- mcp(genstream, alpha=0.05, delta="adaptive", cp=0.99)

#should find confidence interval of length 0.02 centered around 0.05
res
}

