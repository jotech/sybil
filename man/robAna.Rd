\name{robAna}
\alias{robAna}

\encoding{utf8}

\title{
  Robustness Analysis
}

\description{
  Performs robustness analysis for a given model.
}

\usage{
  robAna(model,
         ctrlreact,
         numP = 20,
         lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
         solver = SYBIL_SETTINGS("SOLVER"),
         method = SYBIL_SETTINGS("METHOD"),
         solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM"),
         fld = FALSE, verboseMode = 2, ...)
}

\arguments{
  \item{model}{
    An object of class \code{\link{modelorg}}.
  }
  \item{ctrlreact}{
    An object of class \code{\linkS4class{reactId}}, character or integer.
    Specifies the control reaction -- the parameter to vary.
  }
  \item{numP}{
    The number of points to analyse.
  }
  \item{lpdir}{
    Character value, direction of optimisation.  Can be set to \code{"min"} or
    \code{"max"}.\cr
    Default: \code{SYBIL_SETTINGS("OPT_DIRECTION")}.
  }
  \item{solver}{
    Single character value.  The solver to use.  See
    \code{\link{SYBIL_SETTINGS}} for possible values.\cr
    Default: \code{SYBIL_SETTINGS("SOLVER")}.
  }
  \item{method}{
    Single character value.  The optimization algorithm to use.  Possible values
    depend on the setting in \code{solver}.  See \code{\link{SYBIL_SETTINGS}}
    for possible values.\cr
    Default: \code{SYBIL_SETTINGS("METHOD")}.
  }
  \item{solverParm}{
    A data frame containing parameters for the specified solver.\cr
    Default: SYBIL_SETTINGS("SOLVER_CTRL_PARM").
  }
  \item{fld}{
    Boolean.  Save the resulting flux distribution.\cr
    Default: \code{FALSE}
  }
  \item{verboseMode}{
    An integer value indicating the amount of output to stdout:
    0: nothing, 1: status messages, 2: like 1 plus a progress indicator,
    3: a table containing the control reaction id and the corresponding
    flux values.\cr
    Default: \code{2}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{simpleFBA}}.
  }
}

\details{
  The function \code{robAna} performs a robustness analysis with a given
  model. The flux of \code{ctrlreact} will be varied in \code{numP} steps
  between the maximum and minimum value the flux of \code{ctrlreact} can reach.
  For each of the \code{numP} datapoints the followong lp problem is solved
  \deqn{%
      \begin{array}{rll}%
          \max            & \mbox{\boldmath$c$\unboldmath}^{\mathrm{T}}
                                         \mbox{\boldmath$v$\unboldmath} \\[1ex]
          \mathrm{s.\,t.} & \mbox{\boldmath$Sv$\unboldmath} = 0         \\[1ex]
                          & v_j = c_k                                   \\[1ex]
                          & \alpha_i \leq v_i \leq \beta_i
                          & \quad \forall i \in \{1, \ldots, n\}, i \neq j\\[1ex]
      \end{array}%
  }{
      max  c^T v
      s.t. Sv = 0
           v_j = c_k
           a_i <= v_i <= b_i  for i = 1, ..., n, i != j
  }
  with \eqn{\bold{S}}{S} beeing the stoichiometric matrix, \eqn{\alpha_i}{a_i}
  and \eqn{\beta_i}{b_i} beeing the lower and upper bounds for flux (variable)
  \eqn{i}.  The total number of variables of the optimization problem is denoted
  by \eqn{n}. The parameter \eqn{c_k} is varied \code{numP} times in the range
  of \eqn{v_{j,\mathrm{min}}}{v_i,min} to \eqn{v_{j,\mathrm{max}}}{v_i,max}.
  The result of the optimization is returned as object of class
  \code{\linkS4class{optsol_robAna}} containing the objective value for each
  datapoint.
  
  The extreme points of the range for \code{ctrlreact} are calculated via flux
  balance analysis (see also \code{\link{simpleFBA}}) with the objective
  function being minimization and maximization of the flux through
  \code{ctrlreact}.
}

\value{
  An object of class \code{\linkS4class{optsol_robAna}}.
}

\references{
  Becker, S. A., Feist, A. M., Mo, M. L., Hannum, G., Palsson, B. Ø.  and
  Herrgard, M. J. (2007) Quantitative prediction of cellular metabolism with
  constraint-based models: the COBRA Toolbox. \emph{Nat Protoc} \bold{2},
  727--738.

  Schellenberger, J., Que, R., Fleming, R. M. T., Thiele, I., Orth, J. D.,
  Feist, A. M., Zielinski, D. C., Bordbar, A., Lewis, N. E., Rahmanian, S.,
  Kang, J., Hyduke, D. R. and Palsson, B. Ø. (2011) Quantitative prediction of
  cellular metabolism with constraint-based models: the COBRA Toolbox v2.0.
  \emph{Nat Protoc} \bold{6}, 1290--1307.

  Bernhard Ø. Palsson (2006).
  \emph{Systems Biology: Properties of Reconstructed Networks}.
  Cambridge University Press.
}

\author{Gabriel Gelius-Dietrich}

\examples{
  data(Ec_core)
  rb <- robAna(Ec_core, ctrlreact = "EX_o2")
  plot(rb)
}

\keyword{optimize}

