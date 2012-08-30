\name{getSamples-methods}
\alias{getSamples}
\docType{methods}
\alias{getSamples-methods}
\alias{getSamples,mmctSamplerGeneric-method}
\alias{getSamples,mmctSampler-method}
\title{Methods for Function `cont' in class `mmctestres', Package `simctest'}
\description{
  Function to request further samples from certain hypotheses.
}
\usage{
 getSamples(obj, ind, n)
}
\arguments{
  \item{obj}{object of type "mmctSampler" derived from class "mmctSamplerGeneric".
  }
  \item{ind}{vector containing the indices of hypotheses
for which further samples are requested.
  }
  \item{n}{vector containing number of further samples for each hypothesis in vector `ind'.
  }
}
\section{Methods}{
  \describe{
    \item{getSamples(obj, ind, n)}{ works with object of type "mmctSampler" derived from class "mmctSamplerGeneric". }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n[i])<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  samples <- getSamples(i, c(1,2), c(2,2));
}
\keyword{methods}
