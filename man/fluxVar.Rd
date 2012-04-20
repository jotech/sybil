\name{fluxVar}
\alias{fluxVar}

\encoding{utf8}

\title{
  Flux Variability Analysis
}
\description{
  Performs flux variablilty analysis for a given model.
}
\usage{
  fluxVar(model, react, percentage = 100,
          tol = SYBIL_SETTINGS("TOLERANCE"),
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
  \item{react}{
    An object of class \code{\linkS4class{reactId}}, character or integer.
    Specifies the fluxes (variables) to analyse.\cr
    If left empty: \code{react_id(model)}.
  }
  \item{percentage}{
    Consider solutions with \code{x} percent of the optimal solution.\cr
    Default: \code{100}
  }
  \item{tol}{
    Tolerance value.\cr
    Default: \code{SYBIL_SETTINGS("TOLERANCE")}.
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
    3: a table containing the reaction id's and the corresponding
    min max values.\cr
    Default: \code{2}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{simpleFBA}}.
  }
}

\details{
  The function \code{fluxVar} performs a flux variability analysis with a given
  model. The minimum and maximum flux values for each reaction in the model
  are calculated, which  still support a given optimal functional state
  \eqn{Z_{\mathrm{opt}}}{Z[opt]}. For each flux \eqn{i} two linear programming
  problems are solved
  \deqn{%
      \begin{array}{rll}%
          \max \textrm{ or } \min & v_i                                 \\[1ex]
          \mathrm{s.\,t.} & Z = Z_{\mathrm{opt}}                        \\[1ex]
                          & \mbox{\boldmath$Sv$\unboldmath} = 0         \\[1ex]
                          & \alpha_i \leq v_i \leq \beta_i
                            & \quad \forall i \in \{1, \ldots, n\}      \\[1ex]
      \end{array}%
  }{
      max or min v[i]
      s.t. Z = Z[opt]
           Sv = 0
           a_i <= v_i <= b_i  for i = 1, ..., n
  }
  with \eqn{\bold{S}}{S} beeing the stoichiometric matrix, \eqn{\alpha_i}{a_i}
  and \eqn{\beta_i}{b_i} beeing the lower and upper bounds for flux (variable)
  \eqn{i}.  The total number of variables of the optimization problem is denoted
  by \eqn{n}.  
  The result of the optimization is returned as object of class
  \code{\linkS4class{optsol_fluxVar}} containing the range of each flux still
  supporting the given optimal state.
  
  The optimal state \eqn{Z_{\mathrm{opt}}}{Z[opt]} is calculated via flux
  balance analysis (see also \code{\link{simpleFBA}}). The objective function
  here is the one given in the model.
}

\value{
  An object of class \code{\linkS4class{optsol_fluxVar}}.
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
  fv <- fluxVar(Ec_core)
  plot(fv)
}

\keyword{optimize}

