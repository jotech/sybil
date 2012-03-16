\name{optsol-class}

\docType{class}

\alias{optsol-class}
\alias{optsol}

\alias{checkStat,optsol-method}
\alias{checkStat}

\alias{fluxdist,optsol-method}
\alias{fluxdist}
\alias{fluxdist<-,optsol-method}
\alias{fluxdist<-}
\alias{fluxes,optsol-method}
\alias{fluxes}
\alias{fluxes<-,optsol-method}
\alias{fluxes<-}
\alias{lp_dir,optsol-method}
\alias{lp_dir}
\alias{lp_dir<-,optsol-method}
\alias{lp_dir<-}
\alias{lp_num_cols,optsol-method}
\alias{lp_num_cols}
\alias{lp_num_cols<-,optsol-method}
\alias{lp_num_cols<-}
\alias{lp_num_rows,optsol-method}
\alias{lp_num_rows}
\alias{lp_num_rows<-,optsol-method}
\alias{lp_num_rows<-}
\alias{lp_obj,optsol-method}
\alias{lp_obj}
\alias{lp_obj<-,optsol-method}
\alias{lp_obj<-}
\alias{lp_ok,optsol-method}
\alias{lp_ok}
\alias{lp_ok<-,optsol-method}
\alias{lp_ok<-}
\alias{lp_stat,optsol-method}
\alias{lp_stat}
\alias{lp_stat<-,optsol-method}
\alias{lp_stat<-}
\alias{method,optsol-method}
\alias{method}
\alias{method<-,optsol-method}
\alias{method<-}
\alias{nfluxes,optsol-method}
\alias{nfluxes}
\alias{num_of_prob,optsol-method}
\alias{num_of_prob}
\alias{num_of_prob<-,optsol-method}
\alias{num_of_prob<-}
\alias{obj_function,optsol-method}
\alias{obj_function}
\alias{obj_function<-,optsol-method}
\alias{obj_function<-}
\alias{react_id,optsol-method}
\alias{react_id<-,optsol-method}
\alias{show,optsol-method}
\alias{solver,optsol-method}
\alias{solver}
\alias{solver<-,optsol-method}
\alias{solver<-}
\alias{length,optsol-method}
\alias{hist,optsol-method}

\encoding{utf8}

\title{Class optsol}

\description{
  Structure of the class \code{optsol}. This class is extended by
  other classes and will not be useed as is. Objects of class
  \code{optsol} are used as "base class".
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optsol}:

  \code{test <- optsol(solver)}.
  \describe{
    \item{\code{solver}:}{Object of class \code{"character"}
      indicating the used solver.
    }
  }
}

\section{Slots}{
  \describe{
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

\section{Methods}{
  \describe{
    \item{\code{allGenes<-}:}{\code{signature(object = "optsol")}:
      sets the \code{allGenes} slot.
    }
    \item{\code{allGenes}:}{\code{signature(object = "optsol")}:
      gets the \code{allGenes} slot.
    }
    \item{\code{checkStat}:}{\code{signature(object = "optsol")}:
      returns the indices of problems with a non optimal solution status.
    }
    \item{\code{fluxdist<-}:}{\code{signature(object = "optsol")}:
      sets the \code{fluxdist} slot.
    }
    \item{\code{fluxdist}:}{\code{signature(object = "optsol")}:
      gets the \code{fluxdist} slot.
    }
    \item{\code{fluxes<-}:}{\code{signature(object = "optsol")}:
      sets the \code{fluxes} slot of slot \code{fluxdist}.
    }
    \item{\code{fluxes}:}{\code{signature(object = "optsol")}:
      gets the \code{fluxes} slot of slot \code{fluxdist}.
    }
    \item{\code{lp_dir<-}:}{\code{signature(object = "optsol")}:
      sets the \code{lp_dir} slot.
    }
    \item{\code{lp_dir}:}{\code{signature(object = "optsol")}:
      gets the \code{lp_dir} slot.
    }
    \item{\code{lp_num_cols<-}:}{\code{signature(object = "optsol")}:
      sets the \code{lp_num_cols} slot.
    }
    \item{\code{lp_num_cols}:}{\code{signature(object = "optsol")}:
      gets the \code{lp_num_cols} slot.
    }
    \item{\code{lp_num_rows<-}:}{\code{signature(object = "optsol")}:
      sets the \code{lp_num_rows} slot.
    }
    \item{\code{lp_num_rows}:}{\code{signature(object = "optsol")}:
      gets the \code{lp_num_rows} slot.
    }
    \item{\code{lp_obj<-}:}{\code{signature(object = "optsol")}:
      sets the \code{lp_obj} slot.
    }
    \item{\code{lp_obj}:}{\code{signature(object = "optsol")}:
      gets the \code{lp_obj} slot.
    }
    \item{\code{lp_ok<-}:}{\code{signature(object = "optsol")}:
      sets the \code{lp_ok} slot.
    }
    \item{\code{lp_ok}:}{\code{signature(object = "optsol")}:
      gets the \code{lp_ok} slot.
    }
    \item{\code{lp_stat<-}:}{\code{signature(object = "optsol")}:
      sets the \code{lp_stat} slot.
    }
    \item{\code{lp_stat}:}{\code{signature(object = "optsol")}:
      gets the \code{lp_stat} slot.
    }
    \item{\code{method<-}:}{\code{signature(object = "optsol")}:
      sets the \code{method} slot.
    }
    \item{\code{method}:}{\code{signature(object = "optsol")}:
      gets the \code{method} slot.
    }
    \item{\code{nfluxes}:}{\code{signature(object = "optsol")}:
      gets the number of flux distributions.
    }
    \item{\code{num_of_prob<-}:}{\code{signature(object = "optsol")}:
      sets the \code{num_of_prob} slot.
    }
    \item{\code{num_of_prob}:}{\code{signature(object = "optsol")}:
      gets the \code{num_of_prob} slot.
    }
    \item{\code{obj_function<-}:}{\code{signature(object = "optsol")}:
      sets the \code{obj_function} slot.
    }
    \item{\code{obj_function}:}{\code{signature(object = "optsol")}:
      gets the \code{obj_function} slot.
    }
    \item{\code{react_id<-}:}{\code{signature(object = "optsol")}:
      sets the \code{react_id} slot.
    }
    \item{\code{react_id}:}{\code{signature(object = "optsol")}:
      gets the \code{react_id} slot.
    }
    \item{\code{show}:}{\code{signature(object = "optsol")}:
      prints a summary of the content of instance of class \code{optsol}.
    }
    \item{\code{solver<-}:}{\code{signature(object = "optsol")}:
      sets the \code{solver} slot.
    }
    \item{\code{solver}:}{\code{signature(object = "optsol")}:
      gets the \code{solver} slot.
    }
    \item{\code{length}:}{\code{signature(object = "optsol")}:
      returns the number of optimizations.
    }
    \item{\code{hist}:}{
      \code{signature(x = "optsol")}: plots a histogram of the
      values of the objective function in optimal state.
    }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}},
  \code{\linkS4class{optsol_simpleFBA}}
  \code{\linkS4class{optsol_fluxdel}},
  \code{\linkS4class{optsol_genedel}},
  \code{\linkS4class{optsol_doublegenedel}},
  \code{\linkS4class{optsol_doublefluxdel}},
  \code{\linkS4class{optsol_blockedReact}},
  \code{\linkS4class{optsol_robAna}} and
  \code{\linkS4class{optsol_fluxVar}}
}

\examples{
  showClass("optsol")
}

\keyword{classes}
