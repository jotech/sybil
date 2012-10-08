\name{optsol_optimizeProb-class}

\Rdversion{1.1}
\encoding{utf8}

\docType{class}

\alias{optsol_optimizeProb-class}
\alias{postProc,optsol_optimizeProb-method}
\alias{postProc<-,optsol_optimizeProb-method}
\alias{preProc,optsol_optimizeProb-method}
\alias{preProc<-,optsol_optimizeProb-method}
\alias{preProc}
\alias{preProc<-}
\alias{postProc}
\alias{postProc<-}

\title{Class \code{"optsol_optimizeProb"}}

\description{
  Structure of the class \code{"optsol_optimizeProb"}. Objects of that class are
  returned by the function \code{\link{optimizeProb}} with the argument
  \code{retOptSol} set to \code{TRUE}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("optsol_optimizeProb", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{preProc}:}{Object of class \code{"ppProc"}
      containing the results of pre-processing.
      See also \code{\link{optimizeProb}}.
    }
    \item{\code{postProc}:}{Object of class \code{"ppProc"}
      containing the results of post-processing.
      See also \code{\link{optimizeProb}}.
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
    \item{\code{obj_coef}:}{Object of class \code{"numeric"}
      containing the objective coefficients of the used model.
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
\section{Extends}{
  Class \code{"\linkS4class{optsol}"}, directly.
}

\section{Methods}{
  \describe{
    \item{\code{preProc}:}{\code{signature(object = "optsol_optimizeProb")}:
      gets the \code{preProc} slot.
    }
    \item{\code{preProc<-}:}{\code{signature(object = "optsol_optimizeProb")}:
      sets the \code{preProc} slot.
    }
    \item{\code{postProc}:}{\code{signature(object = "optsol_optimizeProb")}:
      gets the \code{postProc} slot.
    }
    \item{\code{postProc<-}:}{\code{signature(object = "optsol_optimizeProb")}:
      sets the \code{postProc} slot.
    }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}},
  \code{\linkS4class{optsol}},
  \code{\linkS4class{optsol_genedel}} and
  \code{\linkS4class{optsol_fluxdel}}
}


\examples{
  showClass("optsol_optimizeProb")
}

\keyword{classes}
