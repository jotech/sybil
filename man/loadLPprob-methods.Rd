\name{loadLPprob-methods}

\docType{methods}
\encoding{utf8}

\alias{loadLPprob}
\alias{loadLPprob-methods}
\alias{loadLPprob,optObj_clpAPI-method}
\alias{loadLPprob,optObj_cplexAPI-method}
\alias{loadLPprob,optObj_glpkAPI-method}
\alias{loadLPprob,optObj_lpSolveAPI-method}

\title{Load Data to Optimization Problem}

\description{
  Load data to the problem object (extending class
  \code{\linkS4class{optObj}}). Use this method to generate problem
  objects.
}

\usage{
\S4method{loadLPprob}{optObj_clpAPI}(lp,
           nCols, nRows, mat, ub, lb, obj, rlb, rtype,
           lpdir = "max", rub = NULL, ctype = NULL)

\S4method{loadLPprob}{optObj_cplexAPI}(lp,
           nCols, nRows, mat, ub, lb, obj, rlb, rtype,
           lpdir = "max", rub = NULL, ctype = NULL)

\S4method{loadLPprob}{optObj_glpkAPI}(lp,
           nCols, nRows, mat, ub, lb, obj, rlb, rtype,
           lpdir = "max", rub = NULL, ctype = NULL)

\S4method{loadLPprob}{optObj_lpSolveAPI}(lp,
           nCols, nRows, mat, ub, lb, obj, rlb, rtype,
           lpdir = "max", rub = NULL, ctype = NULL)
}

\section{Methods}{
  \describe{
    \item{\code{signature(lp = "optObj_clpAPI")}}{
      method to use with package \pkg{clpAPI}.
    }
    \item{\code{signature(lp = "optObj_cplexAPI")}}{
      method to use with package \pkg{cplexAPI}.
    }
    \item{\code{signature(lp = "optObj_glpkAPI")}}{
      method to use with package \pkg{glpkAPI}.
    }
    \item{\code{signature(lp = "optObj_lpSolveAPI")}}{
      method to use with package \pkg{lpSolveAPI}.
    }
  }
}

\arguments{
  \item{lp}{
    An object of class
    \code{\linkS4class{optObj_clpAPI}},
    \code{\linkS4class{optObj_cplexAPI}},
    \code{\linkS4class{optObj_glpkAPI}} or
    \code{\linkS4class{optObj_lpSolveAPI}}.
  }
  \item{nCols}{
    Number of columns (variables) of the constraint matrix.
  }
  \item{nRows}{
    Number of rows (constraints) of the constraint matrix.
  }
  \item{mat}{
    An object of class \code{\link{Matrix}}. The constraint matrix of the
    problem object. The number of columns in \code{mat} must be \code{nCols}
    and the number of rows in \code{mat} must be \code{nRows}.
  }
  \item{ub}{
    A numeric vector of length \code{nCols} giving the upper bounds of the
    variables of the problem object.
  }
  \item{lb}{
    A numeric vector of length \code{nCols} giving the lower bounds of the
    variables of the problem object.
  }
  \item{obj}{
    A numeric vector of length \code{nCols} giving the objective coefficients
    of the variables of the problem object.
  }
  \item{rlb}{
    A numeric vector of length \code{nRows} giving the right hand side of the
    problem object. If argument \code{rub} is not \code{NULL}, \code{rlb}
    contains the lower bounds of the constraints of the problem object.
  }
  \item{rtype}{
    A character vector of length \code{nRows} giving the constraint type:
    \tabular{lll}{
      \code{"F"}: \tab free variable, no upper or lower bound \tab \eqn{-\infty < x < \infty}{-INF < x < INF}               \cr
      \code{"L"}: \tab variable with lower bound              \tab \eqn{\mathrm{lb} \leq x < \infty}{lb <= x < INF}         \cr
      \code{"U"}: \tab variable with upper bound              \tab \eqn{-\infty < x \leq \mathrm{ub}}{-INF < x <= ub}       \cr
      \code{"D"}: \tab double-bounded variable                \tab \eqn{\mathrm{lb} \leq x \leq \mathrm{ub}}{lb <= x <= ub} \cr
      \code{"E"}: \tab fixed variable                         \tab \eqn{\mathrm{lb} = x = \mathrm{ub}}{lb = x = ub}         \cr
      \code{"R"}: \tab ranged constraint                                                                                    \cr
    }
    Value \code{"R"} is not available for every solver software. Check
    documentation of the solver software if ranged constraints are supported.
    If \code{rtype[i]} is not \code{"F"}, \code{"L"}, \code{"U"}, \code{"D"},
    \code{"E"} or \code{"R"}, the value of \code{rtype[i]} will be set to
    \code{"E"}.
  }
  \item{lpdir}{
    Single character string containing the direction of optimization.
    Can be set to \code{"min"} or \code{"max"}.\cr
    Default: \code{"max"}.
  }
  \item{rub}{
    A numeric vector of length \code{nRows} giving the right hand side of the
    problem object. If not \code{NULL}, it contains the upper bounds of the
    constraints of the problem object.\cr
    Default: \code{NULL}.
  }
  \item{ctype}{
    A character vector of length \code{nCols} giving the variable type. If set
    to \code{NULL}, no specific variable type is set, which usually means, all
    variables are treated as continuous variables.\cr
    Default: \code{NULL}.
    \tabular{lll}{
      \code{"C"}: \tab continuous variable  \cr
      \code{"B"}: \tab binary variable     \cr
      \code{"I"}: \tab integer variable    \cr
      \code{"S"}: \tab semi-continuous variable  \cr
      \code{"N"}: \tab semi-integer variable    \cr
    }
    Values \code{"S"} and \code{"N"} are not available for every solver
    software. Check documentation of the solver software if semi-continuous and
    semi-integer variables are supported.
    If \code{ctype[j]} is not \code{"C"}, \code{"B"}, \code{"I"}, \code{"S"},
    or \code{"N"}, the value of \code{ctype[j]} will be set to \code{"C"}.
  }
}

\details{
  Method \code{loadLPprob} can be used any time after a problem object is
  initialized by \code{\link{initProb}}.
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  Superclass \code{\linkS4class{optObj}} and constructor function
  \code{\link{optObj}}.
}

\keyword{methods}
\keyword{optimize}
