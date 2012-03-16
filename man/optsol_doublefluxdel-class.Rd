\name{optsol_doublefluxdel-class}

\docType{class}

\alias{optsol_doublefluxdel-class}
\alias{optsol_doublefluxdel}

\alias{ddelres,optsol_doublefluxdel-method}
\alias{ddelres}
\alias{delmat,optsol_doublefluxdel-method}
\alias{delmat<-,optsol_doublefluxdel-method}
\alias{delmat}
\alias{delmat<-}
\alias{lethal,optsol_doublefluxdel-method}
\alias{lethal<-,optsol_doublefluxdel-method}
\alias{lethal}
\alias{lethal<-}

\encoding{utf8}

\title{Class \code{"optsol_doublefluxdel"}}

\description{
  Structure of the class \code{"optsol_doublefluxdel"}. Objects of that class
  are returned by the function \code{\link{doubleFluxDel}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optsol_doublefluxdel}:

  \code{test <- optsol_doublefluxdel(solver, nprob, lpdir, nrows, ncols, delrows, delcols, objf, fld)}.
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
    \item{\code{delrows}:}{Object of class \code{"numeric"}
      number of rows in the matrix od slot \code{delmat}.
    }
    \item{\code{delcols}:}{Object of class \code{"numeric"}
      number of columns in the matrix od slot \code{delmat}.
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
    \item{\code{delmat}:}{Object of class \code{"matrix"} with boolean entries:
      rows react1, columns: react2.
    }
    \item{\code{lethal}:}{Object of class \code{"integer"} containing the
      indices of lethal genes/fluxes upon deletion.
    }
    \item{\code{react_id}:}{Object of class \code{"character"}
      containing the reaction id's.
    }
    \item{\code{allGenes}:}{Object of class \code{"character"}
      containing the gene id's.
    }
    \item{\code{chlb}:}{Object of class \code{"numeric"}
      containing the new (changed) values for the columns lower bounds.
    }
    \item{\code{chub}:}{Object of class \code{"numeric"}
      containing the new (changed) values for the columns upper bounds.
    }
    \item{\code{dels}:}{
      Object of class \code{"matrix"} (one column) containing the position of
      the reaction id's of deleted reactions (fluxes) in the slot
      \code{react_id}.
    }
    \item{\code{algorithm}:}{
      Object of class \code{"character"} containing the used algorithm
	  for metabolic network analysis.
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
  Class \code{"\linkS4class{optsol_fluxdel}"}, directly.
  Class \code{"\linkS4class{optsol_simpleFBA}"}, by class "optsol_fluxdel",
  distance 2.
  Class \code{"\linkS4class{optsol}"}, by class "optsol_fluxdel", distance 3.
}

\section{Methods}{
  \describe{
    \item{\code{delmat}}{\code{signature(object = "optsol_doublefluxdel")}: gets the
      \code{delmat} slot.
    }
    \item{\code{delmat<-}}{\code{signature(object = "optsol_doublefluxdel")}: sets the
      \code{delmat} slot.
    }
    \item{\code{lethal<-}:}{\code{signature(object = "optsol_doublefluxdel")}:
      sets the \code{lethal} slot.
    }
    \item{\code{lethal}:}{\code{signature(object = "optsol_doublefluxdel")}:
      gets the \code{lethal} slot.
    }
    \item{\code{ddelres}:}{
      \code{signature(object = "optsol_doublefluxdel")}(\code{slot, react1,
      react2}): gets the solution for specified reaction pairs.
    }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}},
  \code{\linkS4class{optsol}},
  \code{\linkS4class{optsol_genedel}},
  \code{\linkS4class{optsol_doublegenedel}},
  \code{\linkS4class{optsol_fluxdel}} and
  \code{\linkS4class{optsol_simpleFBA}}
}


\examples{
  showClass("optsol_doublefluxdel")
}

\keyword{classes}
