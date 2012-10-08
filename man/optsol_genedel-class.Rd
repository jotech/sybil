\name{optsol_genedel-class}

\docType{class}

\alias{optsol_genedel-class}

\alias{fluxdels,optsol_genedel-method}
\alias{fluxdels}
\alias{fluxdels<-,optsol_genedel-method}
\alias{fluxdels<-}
\alias{hasEffect,optsol_genedel-method}
\alias{hasEffect}
\alias{hasEffect<-,optsol_genedel-method}
\alias{hasEffect<-}
%\alias{ind2id,optsol_genedel-method}
\alias{deleted,optsol_genedel-method}
%\alias{cpsol,optsol_genedel-method}
%\alias{cpsol}

\encoding{utf8}

\title{Class \code{"optsol_genedel"}}

\description{
  Structure of the class \code{"optsol_genedel"}. Objects of that class are
  returned by the function \code{geneDel}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("optsol_genedel", ...)}.
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
    \item{\code{dels}:}{
      Object of class \code{"matrix"} (one column) containing the position of
      the reaction id's of deleted genes in the slot \code{react_id}.
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
  Class \code{"\linkS4class{optsol_fluxdel}"}, directly.
  Class \code{"\linkS4class{optsol_optimizeProb}"},
  by class "optsol_fluxdel", distance 2.
  Class \code{"\linkS4class{optsol}"}, by class "optsol_fluxdel", distance 3.
}

\section{Methods}{
  \describe{
    \item{\code{fluxdels}:}{
      \code{signature(object = "optsol_genedel")}: gets the \code{fluxdels} slot.
    }
    \item{\code{fluxdels<-}:}{
      \code{signature(object = "optsol_genedel")} sets the \code{fluxdels} slot.
    }
    \item{\code{hasEffect}:}{
      \code{signature(object = "optsol_genedel")}: gets the \code{hasEffect} slot.
    }
    \item{\code{hasEffect<-}:}{
      \code{signature(object = "optsol_genedel")}: sets the \code{hasEffect} slot.
    }
%     \item{\code{ind2id}:}{
%       \code{signature(object = "optsol_genedel")}(\code{slotN}): translates the
%       indices in \code{slotN} into corresponding gene or reaction id's. Possible
%       values for \code{slotN} are: \code{"dels"} or \code{"fluxdels"}.
%     }
    \item{\code{deleted}:}{
      \code{signature(object = "optsol_genedel")}(\code{i}): gets the \code{i}th
      element of the \code{dels} slot.
    }
%     \item{\code{cpsol}:}{
%       \code{signature(object = "optsol_genedel")}(\code{i, z, fld)}: copies solution number \code{i}
%       to solution number \code{z} of solution \code{object}. If \code{fld}
%       is set to \code{TRUE}, also the flux distribution will be copied
%       (default: FALSE).
%     }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{checkOptSol}},
  \code{\linkS4class{optsol}},
  \code{\linkS4class{optsol_fluxdel}} and
  \code{\linkS4class{optsol_optimizeProb}}
}

\examples{
  showClass("optsol_genedel")
}

\keyword{classes}

