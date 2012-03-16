\name{doubleFluxDel}
\alias{doubleFluxDel}

\encoding{utf8}

\title{Double Flux Deletion Experiment}

\description{
  Sets every pairwise combination of fluxes in a model for one time to zero.
}

\usage{
  doubleFluxDel(model, react1, react2, allComb = FALSE,
                lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                solver = SYBIL_SETTINGS("SOLVER"),
                method = SYBIL_SETTINGS("METHOD"),
                fld = FALSE, exex = FALSE, ...)
}

\arguments{
  \item{model}{An object of class \code{modelorg}.}
  \item{react1}{An object of class \code{\link{reactId}} or character or integer, fluxes set to zero (upper and lower bound) (default: \code{react_id(model)}).}
  \item{react2}{An object of class \code{\link{reactId}} or character or integer, fluxes set to zero (upper and lower bound) (default: \code{react_id(model)}).}
  \item{allComb}{Logical: If set to \code{TRUE}, all possible pairwise combinations of react1 and react2 were deleted.
    If set to \code{FALSE}, react1 and react2 must have the same length. Every pair c(react1[i], react2[i]) will be deleted (default: FALSE).}
  \item{lpdir}{Direction of optimization (default: max).}
  \item{solver}{The linear problem solver (default: glpk).}
  \item{method}{The method for the linear problem solver (default: simplex).}
  \item{fld}{Logical: If set to TRUE, the flux distribution for
    each optimization step will be stored. This increases the running
    time significantly (default: FALSE).}
  \item{exex}{Logical: exclude exchange reactions (default: FALSE).}
  \item{...}{Further arguments passed to \code{\link{optimizer}}.}
}

\details{
  The function \code{doubleFluxDel} studies the effect of deletion of pairs of fluxes on
  the current objective function (e.g. biomass production).
  
  The first step is the calculation of a reference solution, which is
  the wild type solution. Then n further optimisations are performed (n:
  number of required optimizations) with \eqn{v_i = 0} for \eqn{i = 1,\dots,n}.
}

\value{
  An object of class \code{\link{optsol_doublefluxdel}}.
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\link{modelorg}}, \code{\linkS4class{optsol}}, \code{\linkS4class{optsol_doublefluxdel}}, \code{\link{checkOptSol}}, \code{\link{optimizer}} and  \code{\link{SYBIL_SETTINGS}}. }

\examples{
  data(Ec_core)
  Ec_ofd <- doubleFluxDel(Ec_core)
}

\keyword{optimize}

