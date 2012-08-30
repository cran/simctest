\name{hBH-methods}
\alias{hBH}
\docType{methods}
\alias{hBH-methods}
\title{Method for class `mcmtest', Package `simctest'}
\description{
  Implementation of the multiple testing procedure by Benjamini-Hochberg.
}
\usage{
 hBH(p, threshold)
}
\arguments{
  \item{p}{object of type "numeric".  }
  \item{threshold}{object of type "numeric".  }
}
\section{Methods}{
  \describe{
    \item{hBH(p, threshold)}{ applies the Benjamini-Hochberg procedure to p-values p with given threshold, returns rejected indices }   
  }
}
\examples{
hBH(runif(10),threshold=0.1)
}
\keyword{methods}
