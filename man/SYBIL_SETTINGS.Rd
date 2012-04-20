\name{SYBIL_SETTINGS}
\alias{SYBIL_SETTINGS}
\alias{SOLVER}
\alias{METHOD}
\alias{TOLERANCE}
\alias{MAXIMUM}
\alias{ALGORITHM}
\alias{OPT_DIRECTION}
\alias{PATH_TO_MODEL}
\alias{SOLVER_CTRL_PARM}

\encoding{utf8}

\title{
  Set and Get SyBiL Parameters
}

\description{
  Manage a set of default parameter settings for SyBiL.
}

\usage{
  SYBIL_SETTINGS(parm, value)
}

\arguments{
  \item{parm}{
    A character string giving the name of the parameter to set.
  }
  \item{value}{
    The corresponding value.
  }
}

\details{
  Typical usages are
  \preformatted{
    SYBIL_SETTINGS(parm, value)
    SYBIL_SETTINGS(parm)
    SYBIL_SETTINGS()
  }
  Possible parameters are:
  \describe{
    \item{\code{"SOLVER"}}{
      The default solver for lp problems. Possible values are depend on your
      installed API package.\cr
      \describe{
          \item{\pkg{glpkAPI}:}{
            \code{"glpk"},
          }
          \item{\pkg{cplexAPI}:}{
            \code{"cplex"},
          }
          \item{\pkg{clpAPI}:}{
            \code{"clp"},
          }
          \item{\pkg{lpSolveAPI}:}{
            \code{"lpSolveAPI"}.
          }
      }
      Default: \code{"glpk"}.
    }
    \item{\code{"METHOD"}}{
      The default method to solve lp problems. Possible values are
      \describe{
          \item{\pkg{glpkAPI}:}{
            \code{"simplex"}, \code{"interior"} or \code{"exact"}.
          }
          \item{\pkg{cplexAPI}:}{
            \code{"lpopt"}, \code{"primopt"} \code{"dualopt"}, \code{"baropt"},
            \code{"hybbaropt"}, \code{"hybnetopt"} or \code{"siftopt"}.
          }
          \item{\pkg{clpAPI}:}{
            \code{"general_solve"}, \code{"inidual"} \code{"iniprimal"},
            \code{"inibarrier"}, \code{"inibarriernoc"}, \code{"idiot"},
            \code{"dual"} or \code{"primal"}.
          }
          \item{\pkg{lpSolveAPI}:}{
            \code{"lp_solve"}.
          }
      }
      Default: \code{"simplex"}.\cr
      If the parameter \code{"SOLVER"} is changed, the corrsponding default
      \code{"METHOD"} is the first one mentioned, e.g. for \code{"cplex"},
      it will be \code{"lpopt"}. This change is done automatically when changing
      the solver. It is not possible, to set a not existing \code{"METHOD"}
      for a particular \code{"SOLVER"}, the corresponding default value will
      be used in such a case.
    }
    \item{\code{"TOLERANCE"}}{
      Tolerance value.\cr
      Default: \code{1E-6}.
    }
    \item{\code{"MAXIMUM"}}{
      Absolute maximum value.\cr
      Default: \code{1000}.
    }
    \item{\code{"ALGORITHM"}}{
      Algorithm to use for e.g. FBA. Can be \code{"FBA"}, \code{"linearMOMA"},
      or \code{"linearMOMA_COBRA"}.\cr
      Default: \code{"FBA"}.
    }
    \item{\code{"OPT_DIRECTION"}}{
      Direction of optimization. Can be \code{"max"} or \code{"min"}.\cr
      Default: \code{"max"}.
    }
    \item{\code{"PATH_TO_MODEL"}}{
      Path to a directory to read or write files.\cr
      Default: \code{"."}.
    }
    \item{\code{"SOLVER_CTRL_PARM"}}{
      A \code{data.frame} giving parameters to the optimizer software
      (e.g. GLPK).\cr
      Default: \code{as.data.frame(NA)}.
    }
  }
}

\value{
  If successful, a set of parameters to SyBiL will be returned.
}

\author{
  Gabriel Gelius-Dietrich

  Maintainer: Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>
}


\examples{
  ## show all current parameters
  SYBIL_SETTINGS()
  
  ## show current setting for "SOLVER"
  SYBIL_SETTINGS("SOLVER")

  ## change current solver to glpkAPI
  SYBIL_SETTINGS("SOLVER", "glpk")
\dontrun{
  ## this needs cplexAPI installed
  ## change current solver to cplexAPI
  SYBIL_SETTINGS("SOLVER", "cplex")
}
}

\keyword{ optimize }
