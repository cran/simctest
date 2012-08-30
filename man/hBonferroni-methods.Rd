\name{hBonferroni-methods}
\alias{hBonferroni}
\docType{methods}
\alias{hBonferroni-methods}
\title{Method for class `mcmtest', Package `simctest'}
\description{
  Implementation of independent (Bonferroni) multiple testing.
}
\usage{
 hBonferroni(p, threshold)
}
\arguments{
  \item{p}{object of type "numeric".  }
  \item{threshold}{object of type "numeric".  }
}
\section{Methods}{
  \describe{
    \item{hBonferroni(p, threshold)}{ performs independent multiple testing using the Bonferroni correction at given threshold, returns rejected indices }   
  }
}
\examples{
hBonferroni(runif(10),threshold=0.1)
}
\keyword{methods}
