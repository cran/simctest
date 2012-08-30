\name{hPC-methods}
\alias{hPC}
\docType{methods}
\alias{hPC-methods}
\title{Method for class `mcmtest', Package `simctest'}
\description{
  Implementation of the multiple testing procedure by Pounds&Cheng.
}
\usage{
 hPC(p, threshold)
}
\arguments{
  \item{p}{object of type "numeric".  }
  \item{threshold}{object of type "numeric".  }
}
\section{Methods}{
  \describe{
    \item{hPC(p, threshold)}{ applies the modification by Pounds&Cheng to p-values p with given threshold, returns rejected indices }   
  }
}
\examples{
hPC(runif(10),threshold=0.1)
}
\keyword{methods}
