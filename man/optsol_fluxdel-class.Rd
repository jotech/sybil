\name{optsol_fluxdel-class}

\docType{class}

\alias{optsol_fluxdel-class}
\alias{optsol_fluxdel}

\alias{react_id,optsol_fluxdel-method}
\alias{react_id<-,optsol_fluxdel-method}
\alias{allGenes,optsol_fluxdel-method}
\alias{allGenes<-,optsol_fluxdel-method}
\alias{chlb,optsol_fluxdel-method}
\alias{chlb<-,optsol_fluxdel-method}
\alias{chlb}
\alias{chlb<-}
\alias{chub,optsol_fluxdel-method}
\alias{chub<-,optsol_fluxdel-method}
\alias{chub}
\alias{chub<-}
\alias{dels,optsol_fluxdel-method}
\alias{dels<-,optsol_fluxdel-method}
\alias{dels}
\alias{dels<-}
\alias{algorithm,optsol_fluxdel-method}
\alias{algorithm<-,optsol_fluxdel-method}
\alias{algorithm}
\alias{algorithm<-}
\alias{ind2id,optsol_fluxdel-method}
\alias{ind2id}
\alias{deleted,optsol_fluxdel-method}
\alias{deleted}
\alias{[,optsol_fluxdel,ANY,ANY,ANY-method}

\encoding{utf8}

\title{Class \code{"optsol_fluxdel"}}

\description{
  Structure of the class \code{"optsol_fluxdel"}. Objects of that class are
  returned by the function \code{\link{oneFluxDel}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optsol_fluxdel}:

  \code{test <- optsol_fluxdel(solver, nprob, lpdir, ncols, nrows, objf, fld, comb = 1)}.
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
    \item{\code{comb}:}{Object of class \code{"numeric"}
      number of combinations to use for flux deletions. This will be the number
      of columns in slot \code{dels}. Default: \code{1}.
    }
  }
}

\section{Slots}{
  \describe{
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
  Class \code{"\linkS4class{optsol_simpleFBA}"}, directly.
  Class \code{"\linkS4class{optsol}"}, by class "optsol_simpleFBA", distance 2.
}

\section{Methods}{
  \describe{
    \item{\code{react_id}:}{
      \code{signature(object = "optsol_fluxdel")}: gets the \code{react_id} slot.
    }
    \item{\code{react_id<-}:}{
      \code{signature(object = "optsol_fluxdel")} sets the \code{react_id} slot.
    }
    \item{\code{allGenes}:}{
      \code{signature(object = "optsol_fluxdel")}: gets the \code{allGenes} slot.
    }
    \item{\code{allGenes<-}:}{
      \code{signature(object = "optsol_fluxdel")} sets the \code{allGenes} slot.
    }
    \item{\code{chlb}:}{
      \code{signature(object = "optsol_fluxdel")}: gets the \code{chlb} slot.
    }
    \item{\code{chlb<-}:}{
      \code{signature(object = "optsol_fluxdel")} sets the \code{chlb} slot.
    }
    \item{\code{chub}:}{
      \code{signature(object = "optsol_fluxdel")}: gets the \code{chub} slot.
    }
    \item{\code{chub<-}:}{
      \code{signature(object = "optsol_fluxdel")}: sets the \code{chub} slot.
    }
    \item{\code{dels}:}{
      \code{signature(object = "optsol_fluxdel")}: gets the \code{dels} slot.
    }
    \item{\code{dels<-}:}{
      \code{signature(object = "optsol_fluxdel")} sets the \code{dels} slot.
    }
    \item{\code{algorithm}:}{
      \code{signature(object = "optsol_fluxdel")}: gets the \code{algorithm} slot.
    }
    \item{\code{algorithm<-}:}{
      \code{signature(object = "optsol_fluxdel")} sets the \code{algorithm} slot.
    }
    \item{\code{ind2id}:}{
      \code{signature(object = "optsol_fluxdel")}(\code{slotN}): translates the
      indices in \code{slotN} into corresponding gene or reaction id's. Possible
      values for \code{slotN} are: \code{"dels"}.
    }
    \item{\code{deleted}:}{
      \code{signature(object = "optsol_fluxdel")}(\code{i}): gets the \code{i}th
      element of the \code{dels} slot.
    }
    \item{\code{[}:}{
      \code{signature(x = "optsol_fluxdel")}: access like a vector.
      \code{x[i]} returns a new object of class \code{optsol_fluxdel} containing
      the \code{i}th deletion experiment.
    }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}},
  \code{\linkS4class{optsol}},
  \code{\linkS4class{optsol_genedel}},
  \code{\linkS4class{optsol_doublegenedel}},
  \code{\linkS4class{optsol_doublefluxdel}} and
  \code{\linkS4class{optsol_simpleFBA}}
}
      
\examples{
  showClass("optsol_fluxdel")
}

\keyword{classes}
