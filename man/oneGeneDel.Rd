\name{oneGeneDel}
\alias{oneGeneDel}

\encoding{utf8}

\title{One Gene Deletion Experiment}

\description{
  Compute in sillico simulations of single-gene knock out mutants.
}

\usage{
  oneGeneDel(model, geneList, ...)
}

\arguments{
  \item{model}{
    An object of class \code{modelorg}.
  }
  \item{geneList}{
    A character vector containing the set of genes to be deleted.  If the
    argument is missing, \code{geneList} is set to \code{allGenes(model)}.
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
  the wild type solution. Then n further optimizations are performed (n:
  number of fluxes) with \eqn{v_i = 0} for \eqn{i = 1,\dots,n}.
}

\value{
  An object of class \code{\link[=optsol_genedel-class]{optsol_genedel}}.
}


\seealso{
  \code{\link{modelorg}}, \code{\linkS4class{optsol}},
  \code{\link[=optsol_genedel-class]{optsol_genedel}},
  \code{\link{checkOptSol}}, \code{\link{optimizer}} and
  \code{\link{SYBIL_SETTINGS}}.
}

\examples{
  data(Ec_core)
  Ec_ofd <- oneFluxDel(Ec_core)
}

\keyword{optimize}

