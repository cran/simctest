\name{getNumber-methods}
\alias{getNumber}
\docType{methods}
\alias{getNumber-methods}
\alias{getNumber,mmctSamplerGeneric-method}
\alias{getNumber,mmctSampler-method}
\title{Methods for Function `cont' in class `mmctestres', Package `simctest'}
\description{
  Function to request number of hypotheses.
}
\usage{
 getNumber(obj)
}
\arguments{
  \item{obj}{object of type "mmctSampler" derived from class "mmctSamplerGeneric".
  }
}
\section{Methods}{
  \describe{
    \item{getNumber(obj)}{ works with object of type "mmctSampler" derived from class "mmctSamplerGeneric". }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n)<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  number <- getNumber(i);
}
\keyword{methods}
