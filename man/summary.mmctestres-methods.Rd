\name{summary.mmctestres-methods}
\alias{summary.mmctestres}
\docType{methods}
\alias{summary.mmctestres-methods}
\alias{summary.mmctestres,mmctestres-method}
\alias{summary.mmctestres,mmctest-method}
\title{Methods for class `mmctestres' and `mmctest', Package `simctest'}
\description{
  Function which shows current estimates of p-values.
}
\usage{
  \method{summary}{mmctestres}(object,...)
}
\arguments{
  \item{object}{object of type `mmctestres' or `mmctest'.
  }
  \item{...}{No further arguments needed. Listed only for compatibility with generic `summary' method.}
}
\section{Methods}{
  \describe{
    \item{summary.mmctestres(object)}{ works with object of type mmctestres or mmctest. }
  }
}
\examples{
  fun <- function(ind,n,data) sapply(1:length(ind), function(i) sum(runif(n)<=data[ind[i]]));
  i <- mmctSampler(fun,num=500,data=runif(500));
  a <- mmctest(h=hBH);
  a <- run(a, i, maxsteps=list(maxnum=1000000,undecided=10));
  summary.mmctestres(a);
}
\keyword{methods}
