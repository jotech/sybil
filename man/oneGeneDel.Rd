\name{oneGeneDel}
\alias{oneGeneDel}

\encoding{utf8}

\title{One Gene Deletion Experiment}

\description{
  Sets every set of fluxes corresponding to a gene for one time to zero.
}

\usage{
  oneGeneDel(model, geneList,
             lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
             solver = SYBIL_SETTINGS("SOLVER"),
             method = SYBIL_SETTINGS("METHOD"), fld = FALSE, ...)
}

\arguments{
  \item{model}{
    An object of class \code{modelorg}.
  }
  \item{geneList}{
    A character vector containing the set of genes to be deleted.  If the
    argument is missing, \code{geneList} is set to \code{allGenes(model)}.
  }
  \item{lpdir}{
    A single character value giving the direction of optimization, can be
    \code{"max"} for maximization or \code{"min"} for minimization.\cr
    Default: \code{lpdir = SYBIL_SETTINGS("OPT_DIRECTION")}.
  }
  \item{solver}{
    A single character value mentioning the desired optimization software,
    can be \code{"glpk"}, \code{"cplex"} or \code{"clp"}.\cr 
    Default: \code{lpdir = SYBIL_SETTINGS("SOLVER")}.
  }
  \item{method}{
    A single character value giving the algorithm the optimization software
    mentioned in \code{solver} should use.\cr
    Default: \code{lpdir = SYBIL_SETTINGS("METHOD")}.
  }
  \item{fld}{
    Logical: If set to TRUE, the flux distribution for each optimization step
    will be stored in the solution object.
    Default: \code{FALSE}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{optimizer}}.  Important ones are
    \code{alg} or \code{solverParm} in order to set parameter values for the
    optimization software.
  }
}

\details{
  The function \code{ondFluxDel} studies the effect of deletion of a
  single gene on the current objective function (e.g. biomass
  production).
  
  The first step is the calculation of a reference solution, which is
  the wild type solution. Then n further optimisations are performed (n:
  number of fluxes) with \eqn{v_i = 0} for \eqn{i = 1,\dots,n}.
}

\value{
  An object of class \code{\link{optsol_genedel}}.
}


\seealso{\code{\link{modelorg}}, \code{\linkS4class{optsol}},
         \code{\linkS4class{optsol_fluxdel}}, \code{\link{checkOptSol}},
         \code{\link{optimizer}} and  \code{\link{SYBIL_SETTINGS}}.
}

\examples{
  data(Ec_core)
  Ec_ofd <- oneFluxDel(Ec_core)
}

\keyword{optimize}

