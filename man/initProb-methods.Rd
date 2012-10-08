\name{initProb-methods}

\docType{methods}
\encoding{utf8}

\alias{initProb}
\alias{initProb-methods}

\alias{initProb,optObj_clpAPI-method}
\alias{initProb,optObj_cplexAPI-method}
\alias{initProb,optObj_glpkAPI-method}
\alias{initProb,optObj_lpSolveAPI-method}

\title{Initialize Problem Object}

\description{
  Initialize Problem Object.
}

\usage{
\S4method{initProb}{optObj_clpAPI}(lp, ...)

\S4method{initProb}{optObj_cplexAPI}(lp, ...)

\S4method{initProb}{optObj_glpkAPI}(lp, ...)

\S4method{initProb}{optObj_lpSolveAPI}(lp, nrows, ncols)
}

\section{Methods}{
  \describe{
    \item{\code{signature(lp = "optObj_clpAPI")}}{
      method to use with package \pkg{optObj_clpAPI}.
    }
    \item{\code{signature(lp = "optObj_cplexAPI")}}{
      method to use with package \pkg{optObj_cplexAPI}.
    }
    \item{\code{signature(lp = "optObj_glpkAPI")}}{
      method to use with package \pkg{optObj_glpkAPI}.
    }
    \item{\code{signature(lp = "optObj_lpSolveAPI")}}{
      method to use with package \pkg{optObj_lpSolveAPI}.
    }
  }
}

\arguments{
  \item{lp}{
    An object extending class \code{\link[=optObj-class]{optObj}}.
  }
  \item{nrows}{
    Number of rows (constraints) of the new problem object.
  }
  \item{ncols}{
    Number of columns (variables) of the new problem object.
  }
  \item{...}{
    Further arguments passed to the initialization function of the solver
    package.
  }
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  Superclass \code{\link[=optObj-class]{optObj}} and constructor function
  \code{\link{optObj}}.
}

\keyword{methods}
\keyword{optimize}
