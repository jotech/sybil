\name{optsol_doublegenedel-class}

\docType{class}

\alias{optsol_doublegenedel-class}
\alias{optsol_doublegenedel}

\alias{ind2id,optsol_doublegenedel-method}

\encoding{utf8}

\title{Class \code{"optsol_doublegenedel"}}

\description{
  Structure of the class \code{"optsol_doublegenedel"}. Objects of that class
  are returned by the function \code{\link{doubleGeneDel}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optsol_doublegenedel}:

  \code{test <- optsol_doublegenedel(solver, nprob, lpdir, nrows, ncols, delrows, delcols, objf, fld)}.
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
    \item{\code{fluxdels}:}{
      Object of class \code{"list"} containing the positions of reaction id's
      in the slot \code{react_id} of deleted reactions (fluxes).
      \code{fluxdels(optsol_genedel)[[i]][j] = l}: The deletion of gene \code{i}
      requires the deletion of a set of fluxes \code{1..k} (j \eqn{\le} k),
      \code{j} beeing the j'th reaction of that set.
      The value \code{l} points to the corresponding reaction id in the slot
      \code{react_id}.
    }
    \item{\code{hasEffect}:}{
      Object of class \code{"logical"} indicating weather deletion of
      gene \code{i} has an effect or not.
    }
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
  Class \code{"\linkS4class{optsol_genedel}"}, directly.
  Class \code{"\linkS4class{optsol_doublefluxdel}"}, directly.
  Class \code{"\linkS4class{optsol_fluxdel}"}, by class "optsol_genedel",
  distance 2.
  Class \code{"\linkS4class{optsol_simpleFBA}"}, by class "optsol_genedel",
  distance 3.
  Class \code{"\linkS4class{optsol}"}, by class "optsol_genedel", distance 4.
}

\section{Methods}{
  \describe{
    \item{\code{ind2id}:}{
      \code{signature(object = "optsol_doublegenedel")}(\code{slotN}): translates the
      indices in \code{slotN} into corresponding gene or reaction id's. Possible
      values for \code{slotN} are: \code{"dels"}, \code{"lethal"} or \code{"fluxdels"}.
    }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}},
  \code{\linkS4class{optsol}},
  \code{\linkS4class{optsol_genedel}},
  \code{\linkS4class{optsol_doublefluxdel}},
  \code{\linkS4class{optsol_fluxdel}},
  \code{\linkS4class{optsol_simpleFBA}}
}

\examples{
  showClass("optsol_doublegenedel")
}

\keyword{classes}
