\name{optObj-class}

\docType{class}

\alias{optObj-class}
\alias{optObj}
\alias{lpExtPtr-class}
\alias{pointerToProb-class}
\alias{pointerToProb}

\alias{clpPtr-class}
\alias{cplexPtr-class}
\alias{glpkPtr-class}
\alias{cplexPointer-class}

\alias{probType}
\alias{sensitivityAnalysis}

\alias{wrong_type_msg}
\alias{wrong_solver_msg}
\alias{checkSolStat}
\alias{getMeanReturn}
\alias{getMeanStatus}

\alias{addCols,optObj-method}
\alias{addRows,optObj-method}
\alias{addRowsCols,optObj-method}
\alias{addRowsToProb,optObj-method}
\alias{backupProb,optObj-method}
\alias{changeColsBnds,optObj-method}
\alias{changeColsBndsObjCoefs,optObj-method}
\alias{changeObjCoefs,optObj-method}
\alias{changeRowsBnds,optObj-method}
\alias{delProb,optObj-method}
\alias{dim,optObj-method}
\alias{getColPrim,optObj-method}
\alias{getColsLowBnds,optObj-method}
\alias{getColsUppBnds,optObj-method}
\alias{getFluxDist,optObj-method}
\alias{getNumCols,optObj-method}
\alias{getNumNnz,optObj-method}
\alias{getNumRows,optObj-method}
\alias{getObjCoefs,optObj-method}
\alias{getObjDir,optObj-method}
\alias{getObjVal,optObj-method}
\alias{getRedCosts,optObj-method}
\alias{getRowsLowBnds,optObj-method}
\alias{getRowsUppBnds,optObj-method}
\alias{getSolStat,optObj-method}
\alias{getSolverParm,optObj-method}
\alias{initProb,optObj-method}
\alias{loadMatrix,optObj-method}
\alias{loadProblemData,optObj-method}
\alias{loadProblemDataLM,optObj-method}
\alias{loadProblemDataMTF,optObj-method}
\alias{method,optObj-method}
\alias{probType,optObj-method}
\alias{scaleProb,optObj-method}
\alias{sensitivityAnalysis,optObj-method}
\alias{setObjDir,optObj-method}
\alias{setRhsZero,optObj-method}
\alias{setSolverParm,optObj-method}
\alias{solveLp,optObj-method}
\alias{solver,optObj-method}
\alias{writeProb,optObj-method}

\alias{addCols}
\alias{addRows}
\alias{addRowsCols}
\alias{addRowsToProb}
\alias{backupProb}
\alias{changeColsBnds}
\alias{changeColsBndsObjCoefs}
\alias{changeObjCoefs}
\alias{changeRowsBnds}
\alias{delProb}
\alias{getColPrim}
\alias{getColsLowBnds}
\alias{getColsUppBnds}
\alias{getFluxDist}
\alias{getNumCols}
\alias{getNumNnz}
\alias{getNumRows}
\alias{getObjCoefs}
\alias{getObjDir}
\alias{getObjVal}
\alias{getRedCosts}
\alias{getRowsLowBnds}
\alias{getRowsUppBnds}
\alias{getSolStat}
\alias{getSolverParm}
\alias{initProb}
\alias{loadMatrix}
\alias{loadProblemData}
\alias{loadProblemDataLM}
\alias{loadProblemDataMTF}
\alias{scaleProb}
\alias{sensitivityAnalysis}
\alias{setObjDir}
\alias{setRhsZero}
\alias{setSolverParm}
\alias{solveLp}
\alias{writeProb}

\encoding{utf8}

\title{Strukture of Class "optObj"}

\description{
  Structure of the class \code{"optObj"}. Objects of that class are
  returned by the function \code{\link{prepProbObj}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optObj}:

  \code{test <- optObj(solver, method, pType = "lp")}.
  \describe{
    \item{\code{solver}:}{Object of class \code{"character"}
      indicating the solver to use. See \code{\link{SYBIL_SETTINGS}} for
      possible values.
    }
    \item{\code{method}:}{Object of class \code{"character"}
      indicating the method used by \code{solver}. See
      \code{\link{SYBIL_SETTINGS}} for possible values. If missing, the
      default method for \code{solver} is used (see also
      \code{\link{checkDefaultMethod}}).
    }
    \item{\code{pType}:}{Object of class \code{"character"}
      type of optimization problem.\cr
      Default: \code{"lp"} (linear programming). Can be e.g. \code{"mip"} (mixed
      integer programming) or the like.
    }
  }
}

\section{Slots}{
  \describe{
    \item{\code{oobj}:}{Object of class \code{"pointerToProb"} containing
                        a pointer to a problem object (see notes).}
    \item{\code{solver}:}{Object of class \code{"character"} containing
                          the problem solver.}
    \item{\code{method}:}{Object of class \code{"character"} containing
                          the method (algorithm) used by the problem solver.}
    \item{\code{probType}:}{Object of class \code{"character"} containing
                            the problem type.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{addCols}}{\code{signature(lp = "optObj")}
        (\code{lp, j}): adds \code{j} columns to the problem object. }
    \item{\code{addRows}}{\code{signature(lp = "optObj")}(\code{lp, i}):
        adds \code{i} rows to the problem object. }
    \item{\code{addRowsCols}}{\code{signature(lp = "optObj")}(\code{lp, nrows, ncols}):
        adds \code{nrows} rows and \code{ncols} columns to the problem object.
        Do not use this method to delete rows or columns from the problem object.}
    \item{\code{addRowsToProb}}{\code{signature(lp = "optObj")}(\code{lp, i, type, lb, ub, cind, nzval}):
        adds new rows to the problem object. \code{i}: ascending row indices,
        \code{type}: character vector, giving the type of the new constraints:
        \code{F}: free, \code{L}: lower bound, \code{U}: upper bound, \code{D}:
        double bounded, \code{E}: fixed, \code{R}: ranged, \code{lb}: lower bounds,
        \code{ub}: upper bounds, \code{cind}: a list giving the column indices of
        the non zero elements, \code{nzval}: non zero elements.}
    \item{\code{backupProb}}{\code{signature(lp = "optObj")}:
        copies lp to a new problem object. The complete problem (slot
        \code{oobj}) will be copied, not just the pointer. }
    \item{\code{changeColsBnds}}{\code{signature(lp = "optObj")}(\code{lp, j, lb, ub}):
        change column bounds of columns in integer vector \code{j}.
        Do not use this method for model creation, use
        \code{changeColsBndsObjCoefs} instead. \code{j}: containing
        column indices, \code{lb}: lower bounds, \code{up}: upper bounds.
        \code{j}, \code{lb} and \code{ub} must have the same length. }
    \item{\code{changeColsBndsObjCoefs}}{\code{signature(lp = "optObj")}(\code{lp, j, lb, ub, obj_coef}):
        set column bounds and objective coefficients of columns in the integer
        vector \code{j}. Use this method for model creation only.
        \code{j}: containing column indices, \code{lb}: lower bounds,
        \code{ub}: upper bounds, \code{obj_coef}: objective coefficients.
        \code{j}, \code{lb}, \code{ub} and \code{obj_coef} must have the
        same length. }
    \item{\code{changeObjCoefs}}{\code{signature(lp = "optObj")}(\code{lp, j, obj_coef}):
        change the objective coefficients of the columns in the integer vector
        \code{j}. \code{j}: containing column indices, \code{obj_coef}: objective
        coefficients. \code{j} and \code{obj_coef} must have the same length. }
    \item{\code{changeRowsBnds}}{\code{signature(lp = "optObj")}(\code{lp, i, lb, ub}):
        change the row bounds of the rows in the integer vector \code{i}.
        \code{i}: column indices, \code{lb}: lower bounds, \code{up}: upper
        bounds. \code{i}, \code{lb} and \code{ub} must have the same length. }
    \item{\code{delProb}}{\code{signature(lp = "optObj")}:
        deletes (nullifies) the pointer to the problem object. }
    \item{\code{dim}}{\code{signature(x = "optObj")}:
        returns a vector \code{d} of length two with d[1] and d[2] containing
        the number of rows and columns of the constraint matrix.}
    \item{\code{getColPrim}}{\code{signature(lp = "optObj")}(\code{lp, j}):
        gets the primal value of column \code{j} after optimization. }
    \item{\code{getColsLowBnds}}{\code{signature(lp = "optObj")}(\code{lp, j}):
        gets the lower bound of column \code{j}. }
    \item{\code{getColsUppBnds}}{\code{signature(lp = "optObj")}(\code{lp, j}):
        gets the upper bound of column \code{j}. }
    \item{\code{getFluxDist}}{\code{signature(lp = "optObj")}:
        gets the flux distribution after optimization. }
    \item{\code{getNumCols}}{\code{signature(lp = "optObj")}:
        gets the current number of columns in the problem object. }
    \item{\code{getNumNnz}}{\code{signature(lp = "optObj")}:
        gets the current number of non zero elements in the constraint matrix
        of the problem object. }
    \item{\code{getNumRows}}{\code{signature(lp = "optObj")}:
        gets the current number of rows in the problem object. }
    \item{\code{getObjCoefs}}{\code{signature(lp = "optObj")}:
        gets the objective coefficients in the problem object. }
    \item{\code{getObjDir}}{\code{signature(lp = "optObj")}:
        gets the optimization direction. }
    \item{\code{getObjVal}}{\code{signature(lp = "optObj")}:
        gets the value of the objective function after optimization. }
    \item{\code{getRedCosts}}{\code{signature(lp = "optObj")}:
        gets the dual value of all variables after optimization. }
    \item{\code{getRowsLowBnds}}{\code{signature(lp = "optObj")}:
        gets all lower bounds of the problem object. }
    \item{\code{getRowsUppBnds}}{\code{signature(lp = "optObj")}:
        gets all upper bounds of the problem object. }
    \item{\code{getSolStat}}{\code{signature(lp = "optObj")}:
        gets the solution status after optimization. }
    \item{\code{getSolverParm}}{\code{signature(lp = "optObj")}:
        gets the current parameters for the used solver. }
    \item{\code{initProb}}{\code{signature(lp = "optObj")}:
        initializes the problem object. }
    \item{\code{loadMatrix}}{\code{signature(lp = "optObj")}(\code{lp, ...}):
        loads the contraint matrix to the problem object.
        
        Further arguments required for glpk:
        \code{ne}: number of non zero elements, \code{ia}: row indices of the
        non zero elements, \code{ja}: column indices or the non zero elements,
        \code{ar}: the non zero elements. \code{ia}, \code{ja} and \code{ar}
        must be of length \code{ne}.

        Further arguments required for COIN-OR CLP:
        \code{ncols}: number of columns, \code{nrows}: number of rows,
        \code{ia}: nz\$ia, \code{ja}: nz\$ja, \code{ra}: nz\$ra.

        Further arguments required for lpSolveAPI:
        S matrix of the model. }
    \item{\code{loadProblemData}}{\code{signature(lp = "optObj")}(\code{lp, model}):
        loads the problem data to the problem object: sets the constraint
        matrix, the upper and lower bounds for columns, the right hand
        side to zero and the objective coefficients. \code{model}: object of
        class \code{"\linkS4class{optObj}"}}
    \item{\code{loadProblemDataLM}}{\code{signature(lp = "optObj")}
        (\code{lp, model, wtflux, nCols, nRows, COBRAflag = FALSE, lpdir = NA}):
        loads the problem data to the problem object for linear MOMA:
        sets the constraint matrix, the upper and lower bounds for columns, the
        right hand side and the objective coefficients.
        \code{model}: object of class \code{"\linkS4class{modelorg}"}, \code{wtflux}:
        wild type flux distribution, \code{nCols}: number of columns of the
        constraint matrix, \code{nRows}: number of rows of the constraint
        matrix, \code{COBRAflag}: \code{TRUE}, if performing linear MOMA like
        in the COBRA toolbox, \code{lpdir}: optimization direction.
        See also \code{\link{simpleFBA}}}
    \item{\code{loadProblemDataMTF}}{\code{signature(lp = "optObj")}
        (\code{lp, model, wtflux, nCols, nRows}):
        loads the problem data to the problem object for minimization of total
        flux: sets the constraint matrix, the upper and lower bounds for
        columns, the right hand side and the objective coefficients.
        \code{model}: object of class \code{"\linkS4class{modelorg}"}, \code{wtflux}:
        optimal value of the first optimization, \code{nCols}: number of columns
        of the constraint matrix, \code{nRows}: number of rows of the constraint
        matrix. See also \code{\link{simpleFBA}}}
    \item{\code{method}}{\code{signature(lp = "optObj")}:
        gets the \code{method} slot. }
    \item{\code{probType}}{\code{signature(lp = "optObj")}:
        gets the \code{probType} slot. }
    \item{\code{scaleProb}}{\code{signature(lp = "optObj")}(\code{lp, opt}):
        scaling of the constraint matrix. \code{opt}: scaling options depending
        on the specified solver. }
    \item{\code{sensitivityAnalysis}}{\code{signature(lp = "optObj")}:
        perform sensitivity analysis using \code{"\linkS4class{optObj_glpk}"}
        or \code{"\linkS4class{optObj_cplex}"}. }
    \item{\code{setObjDir}}{\code{signature(lp = "optObj")}:
        sets the optimization direction. }
    \item{\code{setRhsZero}}{\code{signature(lp = "optObj")}:
        sets the right hand side of the problem object to zero. }
    \item{\code{setSolverParm}}{\code{signature(lp = "optObj")}:
        sets the current parameters for the used solver. }
    \item{\code{solveLp}}{\code{signature(lp = "optObj")}:
        solves the optimization problem with the solver mentioned in slot \code{solver}
        and with the algorithm mentioned in slot \code{method}. Returns the exit
        status of the solver. See \code{\link{SYBIL_SETTINGS}} for the available
        solvers and corresponding methods to solve linear optimization problems. }
    \item{\code{solver}}{\code{signature(object = "optObj")}:
        gets the \code{solver} slot. }
    \item{\code{writeProb}}{\code{signature(lp = "optObj")}(\code{lp, fmt, ...}):
        writes problem object to file. \code{fmt}: file type of the output,
        \code{...} further arguments passed to the writing routine depending
        on the used solver (argument \code{fname} is a good candidate here:
        character value giving the file name to write to). }
  }
}

\section{Further usefull Functions}{
  \describe{
    \item{\code{checkSolStat}:}{
      \code{checkSolStat(stat, solver = SYBIL_SETTINGS("SOLVER"))}
      \cr
      Returns the indices of problems with a non optimal solution status.
      \describe{
          \item{\code{stat}}{
            Vector of integer values containing the solution status.
          }
          \item{\code{solver}}{
            Character string specifying the used solver
            (see \code{\link{SYBIL_SETTINGS}}).
          }
      }
    }
    \item{\code{getMeanReturn}:}{
      \code{getMeanReturn(code, solver = SYBIL_SETTINGS("SOLVER"))}
      \cr
      translates the return value (\code{code}) of a solver in a human
      readable string.
    }
    \item{\code{getMeanStatus}:}{
      \code{getMeanStatus(code, solver = SYBIL_SETTINGS("SOLVER"), env = NULL)}
      \cr
      translates the soluton status value (\code{code}) of a solver in a
      human readable string. Argument \code{env} is for use with
      IBM ILOG CPLEX.
    }
    \item{\code{wrong_type_msg}:}{
      \code{wrong_type_msg(lp)}
      \cr
      prints a warning message, if slot \code{oobj} from \code{lp} (an instance
      of class \code{optObj}) does not contain a pointer to a valid solver. See
      also \code{\link{SYBIL_SETTINGS}} for possible solvers.
    }
    \item{\code{wrong_solver_msg}:}{
      \code{wrong_solver_msg(lp, method, printOut = TRUE)}
      \cr
      if \code{printOut == TRUE}, it will print a warning message,
      if \code{method} is not available for \code{solver} in \code{lp}.
    }
  }
}

\details{
  The intention of class \code{optObj} is, to provide a flexible
  user interface to several optimization software products. The
  methods here working on the slot \code{oobj} are interface functions
  to low level functions invoking corresponding C functions.
  Basically, the user has not to care about the nature of the solver,
  or solver-specific functions. That is done by the class.
}

\note{
  The class \code{pointerToProb} contains an external pointer to a problem
  object (usually a C/C++ pointer). This is for \pkg{glpkAPI} an object of class
  \code{\link[glpkAPI]{glpkPtr}}, for \pkg{clpAPI} an object of class
  \code{\link[=externalptr-class]{externalptr}}, for \pkg{lpSolveAPI} an object
  of class \code{lpExtPtr} and for \pkg{cplexAPI} an object of class
  \code{cplexPointer}.
  
  The class \code{cplexPointer} has two slots \code{env} and \code{lp}, each of
  class \code{cplexPtr}. %\code{\link[cplexAPI]{cplexPtr}}.
  To access for example the environment
  pointer from an object of class \code{optObj}, one can write
  \code{lp@oobj@env}.
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\link{prepProbObj}} and \code{\link{SYBIL_SETTINGS}}}

\examples{
  showClass("optObj")
}

\keyword{classes}


