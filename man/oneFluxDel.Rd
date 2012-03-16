\name{oneFluxDel}
\alias{oneFluxDel}

\encoding{utf8}

\title{One Flux Deletion Experiment}

\description{
  Sets every flux in a model for one time to zero. The following optimizations were performed:

  \tabular{rl}{
  maximize   \tab \eqn{c^T \times v}{c^t * v} \cr
  
  subject to \tab \eqn{Sv = b} \cr
    
             \tab \eqn{lb \le v_j \le ub}{lb <= v_j <= ub}, \eqn{j = 0,\dots,i-1,i+1,\dots,n}{j = 0,...,i-1,i+1,...,n} \cr

             \tab \eqn{v_i = 0}
  }
}

\usage{
  oneFluxDel(model, react, lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
             solver = SYBIL_SETTINGS("SOLVER"),
             method = SYBIL_SETTINGS("METHOD"),
             fld = FALSE, ...)
}

\arguments{
  \item{model}{An object of class \code{modelorg}.}
  \item{react}{An object of class \code{\link{reactId}} or character or integer, fluxes set to zero (upper and lower bound) (default: \code{react_id(model)}).}
  \item{lpdir}{Direction of optimization (default: max).}
  \item{solver}{The linear problem solver (default: glpk).}
  \item{method}{The method for the linear problem solver (default: simplex).}
  \item{fld}{Logical: If set to TRUE, the flux distribution for
    each optimization step will be stored. This increases the running
    time significantly (default: FALSE).}
  \item{...}{Further arguments passed to \code{\link{optimizer}}.}
}

\details{
  The function \code{oneFluxDel} studies the effect of deletion of a
  single flux on the current objective function (e.g. biomass
  production).
  
  The first step is the calculation of a reference solution, which is
  the wild type solution. Then n further optimisations are performed (n:
  number of fluxes) with \eqn{v_i = 0} for \eqn{i = 1,\dots,n}.
}

\value{
  An object of class \code{\link{optsol_fluxdel}}.
}

\author{Gabriel Gelius-Dietrich}

\note{
  need to fill in
}

\seealso{\code{\link{modelorg}}, \code{\linkS4class{optsol}}, \code{\linkS4class{optsol_fluxdel}}, \code{\link{checkOptSol}}, \code{\link{optimizer}} and  \code{\link{SYBIL_SETTINGS}}. }

\examples{
  data(Ec_core)
  Ec_ofd <- oneFluxDel(Ec_core)
}

\keyword{optimize}

