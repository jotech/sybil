\name{simpleFBA}

\alias{simpleFBA}
\alias{simpleFBA-deprecated}

\encoding{utf8}

\title{
  Flux Balance Analysis with a Given Model
}

\description{
  Finds an optimal flux distribution and optimal value for the objective
  function for a model.  
}

\usage{
  simpleFBA(model, react = NA, lb = NA, ub = NA, gene = NA,
            lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
            minTotalFlux = FALSE, minDist = FALSE, wtFluxes = NA,
            solver = SYBIL_SETTINGS("SOLVER"),
            method = SYBIL_SETTINGS("METHOD"),
            fld = FALSE, retOptSol = FALSE,
            prCmd = NA, poCmd = NA, prCil = NA, poCil = NA, ...)
}

\arguments{
  \item{model}{
    An object of class \code{\linkS4class{modelorg}} or of class
    \code{\link[=optObj-class]{optObj}}.
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
  \item{gene}{
    A character or integer vector containing gene id's or indices of gene id's
    in \code{\link{allGenes}(model)}. If arguments \code{lb} and/or \code{ub}
    are additionally used (not \code{NA}), upper and lower bounds will be
    applied to all fluxes on which the deletion of the genes given in
    \code{gene} have an effect. In this case, the first value in \code{lb} and
    \code{ub} is used.
    Default: \code{NA}.
  }
  \item{lpdir}{
    Character value, direction of optimization.  Can be set to \code{"min"} or
    \code{"max"}. This parameter is used only if argument \code{model} is
    of class \code{\linkS4class{modelorg}}.\cr
    Default: \code{SYBIL_SETTINGS("OPT_DIRECTION")}.
  }
  \item{minTotalFlux}{
    Boolean, find a flux distribution with minimal total flux leading to an
    optimal value for the objective function (see Details below).\cr
    Default: \code{FALSE}.
  }
  \item{minDist}{
    Boolean, find a flux distribution with a minimal manhattan distance to a
    given flux (see Details below).\cr
    Default: \code{FALSE}.
  }
  \item{wtFluxes}{
    A single numeric value, or a numeric vector of length equal to the number of
    fluxes (variables) in the model (see Details below).\cr
    Default: \code{NA}.
  }
  \item{solver}{
    Single character value.  The solver to use.  See
    \code{\link{SYBIL_SETTINGS}} for possible values. This parameter is used
    only if argument \code{model} is of class \code{\linkS4class{modelorg}}.\cr
    Default: \code{SYBIL_SETTINGS("SOLVER")}.
  }
  \item{method}{
    Single character value.  The optimization algorithm to use.  Possible values
    depend on the setting in \code{solver}.  See \code{\link{SYBIL_SETTINGS}}
    for possible values. This parameter is used only if argument \code{model} is
    of class \code{\linkS4class{modelorg}}.\cr
    Default: \code{SYBIL_SETTINGS("METHOD")}.
  }
  \item{fld}{
    Boolean.  Save the resulting flux distribution.\cr
    Default: \code{FALSE}.
  }
  \item{retOptSol}{
    Boolean.  Return an object of class
    \code{\link[=optsol_optimizeProb-class]{optsol_optimizeProb}} or just a
    list containing the results.\cr
    Default: \code{FALSE}.
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
    Can be used if \code{simpleFBA} is called several times (like in
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
  \item{\dots}{
    Further arguments passed to \code{\link{prepProbObj}}.  Argument
    \code{solverParm} is a good candidate.
  }
}

\details{
  The function \code{simpleFBA} performs a flux balance analysis (FBA) with a
  given model, which is basically the solution of a linear programming problem
  \deqn{%
      \begin{array}{rll}%
          \max            & \mbox{\boldmath$c$\unboldmath}^{\mathrm{T}}
                                         \mbox{\boldmath$v$\unboldmath} \\[1ex]
          \mathrm{s.\,t.} & \mbox{\boldmath$Sv$\unboldmath} = 0         \\[1ex]
                          & \alpha_i \leq v_i \leq \beta_i
                            & \quad \forall i \in \{1, \ldots, n\}      \\[1ex]
      \end{array}%
  }{
      max  c^T v
      s.t. Sv = 0
           a_i <= v_i <= b_i  for i = 1, ..., n
  }
  with \eqn{\bold{S}}{S} beeing the stoichiometric matrix, \eqn{\alpha_i}{a_i}
  and \eqn{\beta_i}{b_i} beeing the lower and upper bounds for flux (variable)
  \eqn{i}.  The total number of variables of the optimization problem is denoted
  by \eqn{n}.  The solution of the optimization is a flux distribution
  maximizing the objective function
  \eqn{
      \mbox{\boldmath$c$\unboldmath}^{\mathrm{T}}
      \mbox{\boldmath$v$\unboldmath}
  }{c^T v} under the a given environment and the assumption of steady state.
  The result of the optimization is returned as list or object of class
  \code{\link[=optsol_optimizeProb-class]{optsol_optimizeProb}} containing the
  optimal value of the objective function and return and status values returned
  by the solver.
  
  When the argument \code{minTotalFlux} is set to \code{TRUE}, a second
  optimization will be performed in order to minimize the sum of the absolute
  values of the fluxes
  \deqn{%
      \begin{array}{rll}%
          \min            &  \begin{minipage}[b]{2.5em}
                             \[
                               \sum_{i=1}^n v_i
                             \]
                             \end{minipage}                              \\[2em]
          \mathrm{s.\,t.} & \mbox{\boldmath$Sv$\unboldmath} = 0          \\[1ex]
                          & \alpha_i \leq v_i \leq \beta_i
                            & \quad \forall i \in \{1, \ldots, n\}       \\[1ex]
                          & \mbox{\boldmath$c$\unboldmath}_{\mathrm{wt}} \geq
                            \mbox{\boldmath$c$\unboldmath}^{\mathrm{T}}
                            \mbox{\boldmath$v$\unboldmath}_{\mathrm{wt}} \\[1ex]
      \end{array}%
  }{
      max  sum v_i  for i = 1, ..., n
      s.t. Sv = 0
           a_i <= v_i <= b_i  for i = 1, ..., n
           c_wt >= c^T v_wt
  }
  with
  \eqn{
      \mbox{\boldmath$c$\unboldmath}^{\mathrm{T}}
      \mbox{\boldmath$v$\unboldmath}_{\mathrm{wt}}
  }{c^T v_wt}
  beeing the optimized value of the objective function of the first
  optimization (wild type solution).  The argument \code{wtFluxes} can be used
  to set a value for a wild type objective function.  In that case, only one
  optimization will be performed.  If \code{wtFluxes} is \code{NA}, a value for
  the wild type will be calculated via FBA.
  
  When the argument \code{minDist} is set to \code{TRUE}, a flux distribution
  leading to an optimal solution as close as possible to a wild type solution
  will be calculated
  \deqn{%
    \begin{array}{rll}%
      \min            & \begin{minipage}[b]{5em}
                        \[ 
                          \sum_{i,j=1}^n
                          \bigl|v_{j,\mathrm{del}} - v_{i,\mathrm{wt}}\bigr|
                        \]
                        \end{minipage}                                 \\[2em]
      \mathrm{s.\,t.} & \mbox{\boldmath$Sv$\unboldmath}_{\mathrm{del}} = 0
                                                                       \\[1ex]
                        & v_i = v_{i,\mathrm{wt}}
                          & \quad \forall i \in \{1, \ldots, n\}       \\[1ex]
                        & \alpha_j \leq v_{j,\mathrm{del}} \leq \beta_j
                          & \quad \forall j \in \{1, \ldots, n\}       \\[1ex]
    \end{array}%
  }{
     max  sum |v_j,del - v_i,wt|  for i,j = 1, ..., n
     s.t. Sv_del = 0
          v_i = vi_i,wt  for i = 1, ..., n
          a_i <= v_i,del <= b_i  for i = 1, ..., n
  }
  Here,
  \eqn{
    \mbox{\boldmath$v$\unboldmath}^{\mathrm{wt}}
  }{v^wt}
  is the wild type flux distribution. This can be set via the argument
  \code{wtFluxes}.  If \code{wtFluxes} is \code{NA} (the default), the
  wild type flux distribution will be calculated by FBA.  If we need to
  calculate a wild type flux distribution and the argument \code{react} is not
  \code{NA}, \code{react} will only be applied to the second optimization.

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

  The paramters \code{solver}, \code{method} and \code{lpdir} are only used if
  argument \code{model} is of class \code{\linkS4class{modelorg}}.
  If \code{model} is of class \code{\link[=optObj-class]{optObj}}, a previous
  call to \code{\link{prepProbObj}} to generate an object of class
  \code{\link[=optObj-class]{optObj}} is assumed. This object already knows
  about this parameters.
}

\note{
  The function \code{simpleFBA} is deprecated and will be defunct in the
  next release. Use \code{\link{optimizeProb}} instead.
}

\value{
  \code{SimpleFBA} returns either an object of class
  \code{\link[=optsol_optimizeProb-class]{optsol_optimizeProb}} if argument
  \code{retOptSol} is set to \code{TRUE}, or simply a list containing the
  results of the optimization:
  \item{ok}{Return value of the optimizer (e.g. \dQuote{solution process was
            successful} or \dQuote{time limit exceeded}).}
  \item{stat}{Status value of the optimization (e.g. \dQuote{solution is
              optimal} or \dQuote{no feasible solution exists}).}
  \item{fluxes}{The resulting flux distribution.}
  \item{preP}{An object of class \code{\link{ppProc}} if a preprocessing command
              was given.}
  \item{postP}{An object of class \code{\link{ppProc}} if a postprocessing
               command was given.}
}

\references{
  Edwards, J. S., Covert, M and Palsson, B. Ø. (2002) Metabolic modelling of
  microbes: the flux-balance approach. \emph{Environ Microbiol} \bold{4},
  133--140.
  
  Edwards, J. S., Ibarra, R. U. and Palsson, B. Ø. (2001) In silico predictions
  of \emph{Escherichia coli} metabolic capabilities are consistent with
  experimental data. \emph{Nat Biotechnol} \bold{19}, 125--130.
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  \code{\link{modelorg}},
  \code{\link[=optsol_optimizeProb-class]{optsol_optimizeProb}},
  \code{\link{prepProbObj}}, \code{\link{checkReactId}},
  \code{\link{SYBIL_SETTINGS}} and \code{\link{ppProc}}.
}

\examples{
\dontrun{
## The examples here require the package glpkAPI to be
## installed. If that package is not available, you have to set
## the argument 'solver' (the default is: solver = "glpkAPI").

## load the example data set
data(Ec_core)

## run simpleFBA(), Ec_sf will be a list
Ec_sf <- simpleFBA(Ec_core)

## run simpleFBA(), Ec_sf will be an object of
## class optsol_optimizeProb
Ec_sf <- simpleFBA(Ec_core, retOptSol = TRUE)

## do FBA, change the upper and lower bounds for the reactions
## "ATPM" and "PFK". Return also the resulting flux distribution
simpleFBA(Ec_core, react = c("ATPM", "PFK"),
          lb = c(3, -3), ub = c(5, 6), fld = TRUE)
          
## do FBA, perform s sensitivity analysis after optimization
simpleFBA(Ec_core, poCmd = list("sensitivityAnalysis"))

## do FBA, write the problem object to file in lp-format
simpleFBA(Ec_core, poCmd = list(c("writeProb", "LP_PROB",
                                  "'lp'", "'Ec_core.lp'"))
)

## do FBA, use "cplexAPI" as lp solver. Get all upper bounds before
## solving the problem. After solving, perform a sensitivity
## analysis and retrieve the reduced costs
simpleFBA(Ec_core, solver = "cplexAPI",
               prCmd = list(c("getColsUppBnds", "LP_PROB", "1:77")),
               poCmd = list("sensitivityAnalysis",
                            c("getDjCPLEX",
                              "LP_PROB@oobj@env",
                              "LP_PROB@oobj@lp",
                              "0", "react_num(Ec_core)-1"
                            )
                       )
)
}
}

\keyword{optimize}

