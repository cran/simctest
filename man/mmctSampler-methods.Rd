\name{mmctSampler-methods}
\alias{mmctSampler}
\docType{methods}
\alias{mmctSampler-methods}

\title{Methods for class mmctSampler', Package `simctest'}
\description{
  Constructor for class 'mmctSampler'.
}
\usage{
 mmctSampler(f, num, data=NULL)
}
\arguments{
  \item{f}{a function f(ind,n,data) which for every hypothesis ind[i]
in vector "ind" returns n[i] new samples and returns the number of
exceedances, where i=1...length(ind).
The data stored in the data slot of class
"mmctSampler" is also passed on to "f".
  }
  \item{num}{number of hypotheses.
  }
  \item{data}{additional slot for data.
  }
}
\section{Methods}{
  \describe{
    \item{mmctSampler(f, num, data)}{ returns object of type 'mmctSampler'
(derived from class 'mmctSamplerGeneric'). }   
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n[i])<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
}
\keyword{methods}
