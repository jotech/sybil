\name{optimizeProb-methods}

\docType{methods}
\encoding{utf8}

\alias{optimizeProb}
\alias{optimizeProb-methods}
\alias{optimizeProb,modelorg-method}
\alias{optimizeProb,sysBiolAlg-method}
\alias{optimizeProb,sysBiolAlg_room-method}

\title{Optimize Problem Object}

\description{
  The function \code{optimizeProb} performs the optimization of a problem
  Object of class \code{\link[=sysBiolAlg-class]{sysBiolAlg}}.
}

\usage{
\S4method{optimizeProb}{modelorg}(object,
             gene = NA,
             react = NA,
             lb = NA,
             ub = NA,
             retOptSol = FALSE,
             MoreArgs = list(),
             ...)

\S4method{optimizeProb}{sysBiolAlg}(object,
             react = NA,
             lb = NA,
             ub = NA,
             obj_coef = NA,
             lpdir = NA,
             resetChanges = TRUE,
             prCmd = NA,
             poCmd = NA,
             prCil = NA,
             poCil = NA)

\S4method{optimizeProb}{sysBiolAlg_room}(object,
             react = NA,
             lb = NA,
             ub = NA,
             prCmd = NA,
             poCmd = NA,
             prCil = NA,
             poCil = NA)
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "modelorg")}}{
      Translates the object of class \code{\linkS4class{modelorg}} into an
      object of class \code{\link[=sysBiolAlg-class]{sysBiolAlg}} and calls
      \code{optimizeProb} again.
    }
    \item{\code{signature(object = "sysBiolAlg")}}{
      Run optimization with the given problem object.
    }
    \item{\code{signature(object = "sysBiolAlg_room")}}{
      Run optimization with the given problem object and take care on specific
      needs for the ROOM algorithm.
      See class \code{\link[=sysBiolAlg_room-class]{sysBiolAlg_room}} for
      details.
    }
  }
}

\arguments{
  \item{object}{
    An object of class \code{\linkS4class{modelorg}} or
    \code{\link[=sysBiolAlg-class]{sysBiolAlg}}.
  }
  \item{gene}{
    A character or integer vector containing gene id's or indices of gene id's
    in \code{\link{allGenes}(model)}. If arguments \code{lb} and/or \code{ub}
    are additionally used (not \code{NA}), upper and lower bounds will be
    applied to all fluxes on which the deletion of the genes given in
    \code{gene} have an effect. In this case, the first value in \code{lb} and
    \code{ub} is used.
    Default: \code{NA}.
  }
  \item{react}{
    An object of class \code{\linkS4class{reactId}}, character or integer.
    Specifies the fluxes (variables) for which to change the upper and lower
    bound (see also arguments \code{lb} and \code{ub}).\cr
    Default: \code{NA}.
  }
  \item{lb}{
    Numeric vector, must have the same length as \code{react}.  Contains the new
    values for the lower bounds of fluxes (variables) mentioned in \code{react}.
    If set to \code{NA}, lower bounds for variables in \code{react} will be set
    to zero.\cr
    Default: \code{NA}.
  }
  \item{ub}{
    Same functionality as \code{lb}, but for upper bounds.\cr
    Default: \code{NA}.
  }
  \item{obj_coef}{
    Numeric vector, must have the same length as \code{react}.  Contains the new
    values for the objective coefficients of fluxes (variables) mentioned in
    \code{react}.  If set to \code{NA}, objective coefficients for variables in
    \code{react} will not be changed.\cr
    Default: \code{NA}.
  }
  \item{lpdir}{
    Character value, direction of optimization.  Can be set to \code{"min"} or
    \code{"max"}.\cr
    Default: \code{SYBIL_SETTINGS("OPT_DIRECTION")}.
  }
  \item{resetChanges}{
    Boolean value. If set to \code{TRUE}, (default) modifications of the
    problem object will be reset to their original values (e.g. changing
    upper and lower bounds for certain reactions). If set to \code{FALSE},
    modifications will stay in the model.\cr
    Default: \code{TRUE}.
  }
  \item{prCmd}{
    A list of preprocessing commands.  See Details below.\cr
    Default: \code{NA}.
  }
  \item{poCmd}{
    A list of postprocessing commands.  See Details below.\cr
    Default: \code{NA}.
  }
  \item{prCil}{
    Can be used if \code{optimizeProb} is called several times (like in
    \code{\link{optimizer}}). The argument \code{prCil} gets the value of the
    loop variable and passes it to the preprocessing function. There, one
    can access it via the keyword \dQuote{LOOP_VAR}. See also
    \code{\link{optimizer}}.\cr
    Default: \code{NA}.
  }
  \item{poCil}{
    Same as \code{prCil}, but for postprocessing.\cr
    Default: \code{NA}.
  }
  \item{retOptSol}{
    Boolean.  Return an object of class
    \code{\link[=optsol_optimizeProb-class]{optsol_optimizeProb}} or
    just a list containing the results.\cr
    Default: \code{FALSE}.
  }
  \item{MoreArgs}{
    Only for the \code{modelorg-method}: a named list containing further
    arguments passed to \code{optimizeProb}. See Details below.
  }
  \item{\dots}{
    Only for the \code{modelorg-method}: further arguments passed to
    \code{\link{sysBiolAlg}}. See Details below.
  }
}

\details{
  The arguments \code{prCmd} and \code{poCmd} can be used to execute \R commands
  working on the problem object.  All commands in \code{prCmd} are executed
  immediately before solving the problem; all commands in \code{poCmd} are
  executed after the problem has been solved.  In all other aspects, the
  arguments work the same.  The value of \code{prCmd} or \code{poCmd} are
  lists of character vectors (each list element is one command).  Each command
  is a character vector and should be built as follows:
  \itemize{
      \item The first element is the name of the function to call.
      \item All other elements are arguments to the function named in the
            first element.
      \item If any argument is character, enclose it in single
            quotes \code{' '}.
      \item Use the keyword \kbd{LP_PROB} in order to refer to the variable name
            of the problem object (object of class
            \code{\link[=optObj-class]{optObj}}).
      \item If the length of the character vector is one, it is treated as a
            function call with the problem object
            (object of class \code{\link[=optObj-class]{optObj}}) as single
            argument.
  }
  The result will be an object of class \code{\link{ppProc}}.
  A few examples for arguments \code{prCmd} or \code{poCmd} (all arguments must
  be lists, see examples section below):
  \preformatted{
    sensitivityAnalysis
  }
  will be translated to the command
  \preformatted{
    sensitivityAnalysis(LP_PROB)
  }
  with \code{LP_PROB} beeing the placeholder for the variable name of the
  problem object. The vector
  \preformatted{
    c("writeProb", "LP_PROB", "'lp'", "'Ec_core.lp'")
  }
  will be translated to the command
  \preformatted{
    writeProb(LP_PROB, 'lp', 'Ec_core.lp')
  }
  The first element will be the function name and the others the arguments to
  that function. The list of commands
  \preformatted{
    list("sensitivityAnalysis",
         c("getDjCPLEX", "LP_PROB@oobj$env",
           "LP_PROB@oobj$lp", "0", "react_num(Ec_core)-1"
         )
    )
  }
  will be translated to the commands
  \preformatted{
    sensitivityAnalysis(LP_PROB)
    getDjCPLEX(LP_PROB@oobj$env, LP_PROB@oobj$lp,
               0, react_num(Ec_core)-1)
  }
  For more information on the usage of \code{prCmd} and \code{poCmd}, see the
  examples section below.

  The method \code{optimizeProb} for class \code{\linkS4class{modelorg}}
  generates a subclass of class \code{\link[=sysBiolAlg-class]{sysBiolAlg}}
  and calls \code{optimizeProb} for that object again. Argument \code{MoreArgs}
  is used to transport arguments to the second \code{optimizeProb} call.
  Argument \code{\dots} instead is used to transport arguments to the constructor
  function \code{\link{sysBiolAlg}}, for example \code{solver}, \code{method}
  and \code{solverParm}.

  Arguments \code{gene}, \code{react}, \code{lb}, \code{ub} and \code{react}
  cause changes in the problem object (object of class
  \code{\link[=optObj-class]{optObj}}, slot \code{problem} of class
  \code{\link[=sysBiolAlg-class]{sysBiolAlg}}). These changes will be reset
  immediately after optimization.
}

\value{
  Calls to \code{optimizeProb} returns either an object of class
  \code{\link[=optsol_optimizeProb-class]{optsol_optimizeProb}}
  if argument \code{retOptSol} is set to \code{TRUE}, or a list containing the
  results of the optimization:
  \item{ok}{
    Return value of the optimizer (e.g. \dQuote{solution process was
    successful} or \dQuote{time limit exceeded}).
  }
  \item{obj}{
    Value of the objective function after optimization.
  }
  \item{stat}{
    Status value of the optimization (e.g. \dQuote{solution is optimal} or
    \dQuote{no feasible solution exists}).
  }
  \item{fluxes}{
    The resulting flux distribution.
  }
  \item{fldind}{
    Pointers to columns (variables) representing a flux (reaction) in the
    original network. The variable \code{fldind[i]} in the solution object
    represents reaction \code{i} in the original network. 
  }
  \item{preP}{
    An object of class \code{\link{ppProc}} if a preprocessing command
    was given.
  }
  \item{postP}{
    An object of class \code{\link{ppProc}} if a postprocessing command
    was given.
  }
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  \code{\link{modelorg}} and \code{\link[=sysBiolAlg-class]{sysBiolAlg}}.
}

\examples{
\dontrun{
## The examples here require the package glpkAPI to be
## installed. If that package is not available, you have to set
## the argument 'solver' (the default is: solver = "glpk").

## load the example data set
data(Ec_core)

## run optimizeProb(), Ec_sf will be a list
Ec_sf <- optimizeProb(Ec_core)

## run optimizeProb(), Ec_sf will be an object of
## class optsol_optimizeProb
Ec_sf <- optimizeProb(Ec_core, retOptSol = TRUE)

## do FBA, change the upper and lower bounds for the reactions
## "ATPM" and "PFK".
optimizeProb(Ec_core, react = c("ATPM", "PFK"),
             lb = c(3, -3), ub = c(5, 6))
          
## do FBA, perform sensitivity analysis after optimization
optimizeProb(Ec_core, MoreArgs = list(poCmd = list("sensitivityAnalysis")))

## do FBA, write the problem object to file in lp-format
optimizeProb(Ec_core,
             MoreArgs = list(poCmd = list(c("writeProb", "LP_PROB",
                                            "'lp'", "'Ec_core.lp'"))))

## do FBA, use "cplex" as lp solver. Get all upper bounds before
## solving the problem. After solving, perform a sensitivity
## analysis and retrieve the reduced costs
optimizeProb(Ec_core, solver = "cplexAPI", MoreArgs = list(
             prCmd = list(c("getColsUppBnds", "LP_PROB", "1:77")),
             poCmd = list("sensitivityAnalysis",
                          c("getDjCPLEX",
                            "LP_PROB@oobj@env",
                            "LP_PROB@oobj@lp",
                            "0", "react_num(Ec_core)-1"
                          )
                     )
             )
)
}
}

\keyword{methods}
\keyword{optimize}
