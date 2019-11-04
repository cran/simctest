\name{NEWS}
\title{News for Package 'simctest'}

\section{Changes in version 2.6 (2019-11-04)}{
  \itemize{
    \item Adjusted code to use p-value buckets as in the final version
  of hte article "Implementing Monte Carlo Tests with P-value Buckets". 
  }
}

\section{Changes in version 2.5}{
  \itemize{
    \item Added functionality for multiple thresholds (for article 'Implementing Monte Carlo Tests with Multiple Thresholds').
    \item Added main functions 'mctest', 'mctest.RL' and 'mctest.simctest' as well as class 'mtestres'. 
    \item NEWS.RD file added.
  }
}

\section{Changes in version 2.4 (2013-12-13)}{
  \itemize{
    \item Simplified examples in Rd files (and vignettes) to speed-up CRAN self-check.
  }
}

\section{Changes in version 2.3 (2013-06-10)}{
  \itemize{
    \item Implemented functionality for article 'QuickMMCTest--Quick Multiple Monte Carlo Testing'
    \item All changes in 'R/mmctest.R': parameter 'thompson' now available to enable 'QuickMMCTest' allocation of new samples.
  }
}

\section{Changes in version 2.2 (2013-02-20)}{
  \itemize{
    \item Added intersecting of confidence intervals in 'mmctest' (file 'R/mmctest.R')
  }
}

\section{Changes in version 2.1 (2012-09-17)}{
  \itemize{
    \item Implemented functionality for article 'MMCTest -- A safe algorithm for multiple Monte Carlo tests'.
    \item Added classes 'mmctest' and 'mmctestres' (in file 'R/mmctest.R') providing
    main functions 'mmctest' and 'mmctSampler' (plus further functionality).
    \item Added manual Rd files in 'man/' and a new vignette in 'inst/simctest-mmctest-intro.Rnw'.
  }
}

\section{Changes in version 1.99-2 (2012-03-19)}{
  \itemize{
    \item New example in manual file 'man/simctest.Rd'.
  }
}

\section{Changes in version 1.99 (2011-11-08)}{
  \itemize{
    \item Added functions 'mcp' and 'mkdeltamid' (computation of the power of a Monte Carlo test).
    \item Added new class 'mcpres' (in R/mcp.R).
  }
}

\section{Changes in version 1.0-1 (2010-02-07)}{
  \itemize{
    \item Completed missing entries in the 'Rd' file manual (folder 'man').
    \item Added the full JASA bib item in inst/CITATION after publication of the underlying JASA article.
  }
}

\section{Changes in version 1.0 (2009-05-28)}{
  \itemize{
    \item Initial submission to CRAN.
  }
}
