\name{loadQobj-methods}

\docType{methods}
\encoding{utf8}

\alias{loadQobj}
\alias{loadQobj-methods}

\alias{loadQobj,optObj_cplexAPI-method}

\title{Load Quadratic Part of the Objective Function to the Optimization Problem}

\description{
 load quadratic part of the objective function to the optimization problem.
}

\usage{
\S4method{loadQobj}{optObj_cplexAPI}(lp, mat)
}

\section{Methods}{
  \describe{
    \item{\code{signature(lp = "optObj_cplexAPI")}}{
      method to use with package \pkg{optObj_cplexAPI}.
    }
  }
}

\arguments{
  \item{lp}{
    An object extending class \code{\linkS4class{optObj}}.
  }
  \item{mat}{
    An object of class \code{\link{Matrix}} containing the quadratic objective
    Matrix \eqn{Q}.
  }
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  Superclass \code{\linkS4class{optObj}} and constructor function
  \code{\link{optObj}}.
}

\keyword{methods}
\keyword{optimize}
