\name{optsol_robAna-class}

\Rdversion{1.1}
\docType{class}

\alias{optsol_robAna-class}
\alias{optsol_robAna}

\alias{ctrlfl,optsol_robAna-method}
\alias{ctrlfl<-,optsol_robAna-method}
\alias{ctrlfl}
\alias{ctrlfl<-}
\alias{ctrlr,optsol_robAna-method}
\alias{ctrlr<-,optsol_robAna-method}
\alias{ctrlr}
\alias{ctrlr<-}
\alias{plot,optsol_robAna,missing-method}

\title{Class \code{"optsol_robAna"}}

\description{
  Structure of the class \code{"optsol_robAna"}. Objects of that class are
  returned by the function \code{\link{robAna}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optsol_robAna}:

  \code{test <- optsol_robAna(solver, method, nprob, lpdir, ncols, nrows, objf, fld, cr, crf)}.
  \describe{
    \item{\code{solver}:}{Object of class \code{"character"}
      indicating the used solver.
    }
    \item{\code{method}:}{Object of class \code{"character"}
      indicating the algorithm used by \code{solver}.
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
    \item{\code{cr}:}{Object of class \code{"reactId"}
      giving the id of the control flux.
    }
    \item{\code{crf}:}{Object of class \code{"numeric"}
      a vector of flux values for the control reaction.
    }
  }
}

\section{Slots}{
  \describe{
    \item{\code{ctrlr}:}{Object of class \code{"reactId"}
      containing the reaction id of the control reaction.
    }
    \item{\code{ctrlfl}:}{Object of class \code{"numeric"}
      containing the control flux values.
    }
    \item{\code{preProc}:}{Object of class \code{"ppProc"}
      containing the results of pre-processing.
      See also \code{\link{simpleFBA}}.
    }
    \item{\code{postProc}:}{Object of class \code{"ppProc"}
      containing the results of post-processing.
      See also \code{\link{simpleFBA}}.
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
    \item{\code{ctrlfl}:}{\code{signature(object = "optsol_robAna")}:
      gets the \code{ctrlfl} slot.
    }
    \item{\code{ctrlfl<-}:}{\code{signature(object = "optsol_robAna")}:
      sets the \code{ctrlfl} slot.
    }
    \item{\code{ctrlr}:}{\code{signature(object = "optsol_robAna")}:
      gets the \code{ctrlr} slot.
    }
    \item{\code{ctrlr<-}:}{\code{signature(object = "optsol_robAna")}:
      sets the \code{ctrlr} slot.
    }
    \item{\code{plot}}{\code{signature(x = "optsol_robAna", y = "missing")}
      (\code{xlab = "Control Flux", ylab = "Objective Function", type = "b",
             pch = 20, fillColorBg = "grey", fillBg = TRUE, ...}):\cr
      plots the optimal values of the objective function vs. the control flux
      values.
      \describe{
        \item{\code{xlab}}{
          label of x-axis, see also \code{\link{par}}.
        }
        \item{\code{ylab}}{
          label of y-axis, see also \code{\link{par}}.
        }
        \item{\code{type}}{
          plot type, see also \code{\link{par}}.
        }
        \item{\code{pch}}{
          how to plot the points, see also \code{\link{par}}.
        }
        \item{\code{fillColorBg}}{
          color of the area below the curve.
        }
        \item{\code{fillBg}}{
          color the area below the curve.
        }
        \item{\code{...}}{
          further graphical parameters to the \code{\link{points}} function.
        }
      }
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
