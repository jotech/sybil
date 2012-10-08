\name{oneFluxDel}
\alias{oneFluxDel}

\encoding{utf8}

\title{One Flux Deletion Experiment}

\description{
  Single reaction (flux) deletion analysis.
}

\usage{
  oneFluxDel(model, react, ...)
}

\arguments{
  \item{model}{
    An object of class \code{modelorg}.
  }
  \item{react}{
    An object of class \code{\link{reactId}} or character or integer, fluxes set
    to zero (upper and lower bound) (default: \code{react_id(model)}).
  }
  \item{\dots}{Further arguments passed to \code{\link{optimizer}}.}
}

\details{
  The function \code{oneFluxDel} studies the effect of deletion of a
  single flux on the current objective function (e.g. biomass
  production).
  
  The first step is the calculation of a reference solution, which is
  the wild type solution. Then n further optimizations are performed (n:
  number of fluxes) with \eqn{v_i = 0} for \eqn{i = 1,\dots,n}.
}

\value{
  An object of class \code{\link[=optsol_fluxdel-class]{optsol_fluxdel}}.
}

\author{Gabriel Gelius-Dietrich}

\note{
  need to fill in
}

\seealso{
  \code{\link{modelorg}}, \code{\linkS4class{optsol}},
  \code{\link[=optsol_fluxdel-class]{optsol_fluxdel}},
  \code{\link{checkOptSol}}, \code{\link{optimizer}} and
  \code{\link{SYBIL_SETTINGS}}.
}

\examples{
  data(Ec_core)
  Ec_ofd <- oneFluxDel(Ec_core)
}

\keyword{optimize}

