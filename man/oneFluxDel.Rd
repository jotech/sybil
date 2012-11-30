\name{oneFluxDel}
\alias{oneFluxDel}

\encoding{utf8}

\title{One Flux Deletion Experiment}

\description{
  Single reaction (flux) deletion analysis.
}

\usage{
  oneFluxDel(model, react = c(1:react_num(model)), ...)
}

\arguments{
  \item{model}{
    An object of class \code{\linkS4class{modelorg}}.
  }
  \item{react}{
    An object of class \code{\link{reactId}} or character or integer containing
    reaction id's to constrain to zero one by one.\cr
    Default: all reactions present in argument \code{model}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{optimizer}}.  Important ones are
    \code{algorithm} in order to set the algorithm to use or \code{solverParm}
    in order to set parameter values for the optimization software.
  }
}

\details{
  The function \code{oneFluxDel} studies the effect of constraining single
  fluxes to zero flux rates on the phenotype of the metabolic network. The
  function performs \eqn{n} optimizations with \eqn{n} being the number of
  reaction id's given in argument \code{react}. Each optimization corresponds
  to the removal of one reaction.
}

\value{
  An object of class \code{\linkS4class{optsol_fluxdel}}.
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{modelorg}}, \code{\linkS4class{optsol}},
  \code{\linkS4class{optsol_fluxdel}},
  \code{\link{checkOptSol}}, \code{\link{optimizer}} and
  \code{\link{SYBIL_SETTINGS}}.
}

\examples{
  data(Ec_core)
  Ec_ofd <- oneFluxDel(Ec_core)
}

\keyword{optimize}

