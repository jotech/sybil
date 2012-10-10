\name{optsol_simpleFBA-class}

\docType{class}

\alias{optsol_simpleFBA-class}
\alias{optsol_simpleFBA}
\alias{preProc,optsol_simpleFBA-method}
\alias{preProc<-,optsol_simpleFBA-method}
\alias{postProc,optsol_simpleFBA-method}
\alias{postProc<-,optsol_simpleFBA-method}

\encoding{utf8}

\title{Class \code{"optsol_simpleFBA"}}

\description{
  Structure of the class \code{"optsol_simpleFBA"}. Objects of that class are
  returned by the function \code{\link{simpleFBA}} with the argument
  \code{retOptSol} set to \code{TRUE}. The class \code{"optsol_simpleFBA"} is
  deprecated, use class \code{"optsol_optimizeProb"} instead.
}
 
\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optsol_simpleFBA}:

  \code{test <- optsol_simpleFBA(solver, nprob, lpdir, ncols, nrows, objf, fld)}.
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
    \item{\code{objf}:}{Object of class \code{"character"}
      a character string describing the objective function.
    }
    \item{\code{fld}:}{Object of class \code{"logical"}
      logical, if \code{TRUE}, memory for flux distributions will be allocated.
    }
  }
}

\section{Slots}{
  \describe{
    \item{\code{preProc}:}{Object of class \code{"ppProc"}
      containing the results of pre-processing.
      See also \code{\link{simpleFBA}}.
    }
    \item{\code{postProc}:}{Object of class \code{"ppProc"}
      containing the results of post-processing.
      See also \code{\link{simpleFBA}}.
    }
    \item{\code{mod_id}:}{Object of class \code{"character"}
      containing the model id of the used model.
    }
    \item{\code{solver}:}{Object of class \code{"character"}
      indicating the used solver.
    }
    \item{\code{method}:}{Object of class \code{"character"}
      indicating the used method.
    }
    \item{\code{algorithm}:}{Object of class \code{"character"}
      containing the name of the algorithm used for optimizations.
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
      containing the optimal values of the objective function after
      optimization. If no flux distribution is available, slot \code{lp_obj}
      contains the cross-product of the objective coefficients in slot
      \code{obj_coef} and the part of the flux distribution in slot
      \code{fluxdist} containing the values representing fluxes in the entire
      metabolic network (slot \code{fldind}).
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
    \item{\code{obj_coef}:}{Object of class \code{"numeric"}
      containing the objective coefficients of the used model. These are not
      necessarily the objective coefficients of the used algorithm. 
    }
    \item{\code{fldind}:}{Object of class \code{"integer"}
      pointers to columns (variables) representing a flux (reaction) in the
      original network. The variable \code{fldind[i]} in the problem object
      represents reaction \code{i} in the original network. 
    }
    \item{\code{fluxdist}:}{Object of class \code{"fluxDistribution"}
      containing the solutons flux distributions.
    }
  }
}

\section{Methods}{
  \describe{
    \item{\code{preProc}:}{\code{signature(object = "optsol_simpleFBA")}:
      gets the \code{preProc} slot.
    }
    \item{\code{preProc<-}:}{\code{signature(object = "optsol_simpleFBA")}:
      sets the \code{preProc} slot.
    }
    \item{\code{postProc}:}{\code{signature(object = "optsol_simpleFBA")}:
      gets the \code{postProc} slot.
    }
    \item{\code{postProc<-}:}{\code{signature(object = "optsol_simpleFBA")}:
      sets the \code{postProc} slot.
    }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{optsol}"}, directly.
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}},
  \code{\linkS4class{optsol}},
  \code{\linkS4class{optsol_genedel}}, and
  \code{\linkS4class{optsol_fluxdel}}
}

\examples{
  showClass("optsol_simpleFBA")
}

\keyword{classes}
