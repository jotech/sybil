\name{getFluxDist-methods}

\docType{methods}
\encoding{utf8}

\alias{getFluxDist}
\alias{getFluxDist-methods}

\alias{getFluxDist,optObj_clpAPI-method}
\alias{getFluxDist,optObj_cplexAPI-method}
\alias{getFluxDist,optObj_glpkAPI-method}
\alias{getFluxDist,optObj_lpSolveAPI-method}

\title{Retrieve Flux Distribution}

\description{
  Get all primal values of variables after optimization (the resulting flux
  distribution).
}

\usage{
\S4method{getFluxDist}{optObj_clpAPI}(lp)

\S4method{getFluxDist}{optObj_cplexAPI}(lp)

\S4method{getFluxDist}{optObj_glpkAPI}(lp)

\S4method{getFluxDist}{optObj_lpSolveAPI}(lp)
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
}

\value{
  A numeric vector containing all primal values (the flux distribution).
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
