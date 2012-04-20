\name{prepProbObj}
\alias{prepProbObj}

\encoding{utf8}

\title{Creates an Instance of Class "optObj"}

\description{
  The function \code{prepProbObj} creates an instance of class
  \code{\linkS4class{optObj}} with respect to the specified
  solver.
}

\usage{
  prepProbObj(model, nCols, nRows,
              wtflux = NA,
              MOMAflag = FALSE,
              COBRAflag = FALSE,
              minTotalFluxFLAG = FALSE,
              minDistFLAG = FALSE,
              scaling = NA,
              solver = SYBIL_SETTINGS("SOLVER"),
              method = SYBIL_SETTINGS("METHOD"),
              lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
              solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM"))
}

\arguments{
  \item{model}{
    An object of class \code{modelorg}.
  }
  \item{nCols}{
    Number of columns (variables) of the problem object.
  }
  \item{nRows}{
    Number of rows (constraints) of the problem object.
  }
  \item{wtflux}{
    A single numeric value, or a numeric vector of length equal to the number of
    fluxes (variables) in the model (see \code{\link{simpleFBA}}).\cr
    Default: \code{NA}.
  }
  \item{MOMAflag}{
    Boolean, prepare problem object in order to perform minimization of
    metabolic adjustment (MOMA, see \code{\link{optimizer}}).\cr
    Default: \code{FALSE}.
  }
  \item{COBRAflag}{
    Boolean, prepare problem object in order to perform minimization of
    metabolic adjustment as in COBRA (MOMA, see \code{\link{optimizer}}).\cr
    Default: \code{FALSE}.
  }
  \item{minTotalFluxFLAG}{
    Boolean, prepare problem object in order to find a flux distribution with
    minimal total flux leading to an optimal value for the objective function
    (see \code{\link{simpleFBA}}).\cr
    Default: \code{FALSE}.
  }
  \item{minDistFLAG}{
    Boolean, prepare problem object in order to find a flux distribution with a
    minimal manhattan distance to a given flux (see \code{\link{simpleFBA}}).\cr
    Default: \code{FALSE}.
  }
  \item{scaling}{
    Scaling options used to scale the constraint matrix. If set to NA, no
    scaling will be performed (see \code{\linkS4class{optObj}}).\cr
    Default: \code{FALSE}.
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
  \item{lpdir}{
    Character value, direction of optimisation.  Can be set to \code{"min"} or
    \code{"max"}.\cr
    Default: \code{SYBIL_SETTINGS("OPT_DIRECTION")}.
  }
  \item{solverParm}{
    A data frame containing parameters for the specified solver.\cr
    Default: \code{SYBIL_SETTINGS("SOLVER_CTRL_PARM")}, (see Details below).
  }
}

\details{
  The function \code{prepProbObj} creates an instance of class
  \code{\linkS4class{optObj}} with respect to the specified
  solver and desired algorithm. The object type (slot oobj) is
  determined by the argument solver. The properties of the
  problem itself are dependent on the specified algorithm
  (e.g. number of rows and columns of the constraint matrix, etc.).
  
  Parameters can be set as data frame:
  \code{solverParm = data.frame(parm1 = val1, parm2 = val2)}
  with \code{parm1} and \code{parm2} beeing the names of
  two different parameters and \code{val1} and \code{val2} the
  corresponding values. For possible parameters and values see
  the documentation of the used solver package (e.g. \pkg{glpkAPI}).
}

\value{
  An object of class \code{\linkS4class{optObj}}, or FALSE.
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
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\linkS4class{optObj}}, \code{\link{optimizer}},
         \code{\link{simpleFBA}} and \code{\link{SYBIL_SETTINGS}}.
}


\keyword{optimize}


