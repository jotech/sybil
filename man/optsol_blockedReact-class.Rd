\name{optsol_blockedReact-class}

\Rdversion{1.1}
\docType{class}

\alias{optsol_blockedReact-class}
\alias{optsol_blockedReact}

\alias{blocked,optsol_blockedReact-method}
\alias{blocked<-,optsol_blockedReact-method}
\alias{blocked}
\alias{blocked<-}
\alias{react,optsol_blockedReact-method}
\alias{react<-,optsol_blockedReact-method}
\alias{react}
\alias{react<-}
\alias{maxSol,optsol_blockedReact-method}
\alias{minSol,optsol_blockedReact-method}
\alias{maxSol}
\alias{minSol}

\title{Class \code{"optsol_blockedReact"}}

\description{
  Structure of the class \code{"optsol_blockedReact"}. Objects of that class are
  returned by the function \code{\link{blockedReact}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optsol_blockedReact}:

  \code{test <- optsol_blockedReact(solver, nprob, lpdir, ncols, nrows, fld)}.
  \describe{
    \item{\code{solver}:}{Object of class \code{"character"}
      indicating the used solver.
    }
    \item{\code{nprob}:}{Object of class \code{"numeric"}
      number of problems to solve.
    }
    \item{\code{lpdir}:}{Object of class \code{"character"}
      direction of optimization.
    }
    \item{\code{ncols}:}{Object of class \code{"numeric"}
      number of variables in the model.
    }
    \item{\code{nrows}:}{Object of class \code{"numeric"}
      number of constraints in the model.
    }
    \item{\code{fld}:}{Object of class \code{"logical"}
      logical, if \code{TRUE}, memory for flux distributions will be allocated.
    }
  }
}

\section{Slots}{
  \describe{
    \item{\code{blocked}:}{Object of class \code{"logical"}
      indicating if a reaction is blocked, or not.
    }
    \item{\code{react}:}{Object of class \code{"reactId"}
      containing the reaction id's.
    }
    \item{\code{solver}:}{Object of class \code{"character"}
      indicating the used solver.
    }
    \item{\code{method}:}{Object of class \code{"character"}
      indicating the used method.
    }
    \item{\code{num_of_prob}:}{Object of class \code{"integer"}
      indicating the number of optimization problems.
    }
    \item{\code{lp_num_cols}:}{Object of class \code{"integer"}
      indicating the number of columns.
    }
    \item{\code{lp_num_rows}:}{Object of class \code{"integer"}
      indicating the number of rows.
    }
    \item{\code{lp_obj}:}{Object of class \code{"numeric"}
      containing the values of the objective function.
    }
    \item{\code{lp_ok}:}{Object of class \code{"integer"}
      containing the exit code of the optimization.
    }
    \item{\code{lp_stat}:}{Object of class \code{"integer"}
      containing the solution status of the optimization.
    }
    \item{\code{lp_dir}:}{Object of class \code{"character"}
      indicating the direction of optimization.
    }
    \item{\code{obj_function}:}{Object of class \code{"character"}
      containing the name of the objective function.
    }
    \item{\code{fluxdist}:}{Object of class \code{"fluxDistribution"}
      containing the solutons flux distributions.
    }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{optsol}"}, directly.
}

\section{Methods}{
  \describe{
    \item{\code{blocked}:}{
      \code{signature(object = "optsol_blockedReact")}: gets the \code{blocked} slot.
    }
    \item{\code{blocked<-}:}{
      \code{signature(object = "optsol_blockedReact")} sets the \code{blocked} slot.
    }
    \item{\code{react}:}{
      \code{signature(object = "optsol_blockedReact")}: gets the \code{react} slot.
    }
    \item{\code{react<-}:}{
      \code{signature(object = "optsol_blockedReact")} sets the \code{react} slot.
    }
    \item{\code{maxSol}:}{
      \code{signature(object = "optsol_blockedReact")}(\code{slot}): returns the
      values in the slot given in \code{slot} for optimizations in \dQuote{max}
      direction.
    }
    \item{\code{minSol}:}{
      \code{signature(object = "optsol_blockedReact")}(\code{slot}): returns the
      values in the slot given in \code{slot} for optimizations in \dQuote{min}
      direction.
    }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}} and \code{\linkS4class{optsol}}
}
      
\examples{
  showClass("optsol_blockedReact")
}

\keyword{classes}
