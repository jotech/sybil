\name{optimizer}
\alias{optimizer}

\encoding{utf8}

\title{Performs Optimizations for the Gene/Flux Deletion Experiments}

\description{
  This function is used by the gene/flux deletion experiments
  (the functions \code{\link{oneFluxDel}}, \code{\link{oneGeneDel}},
  \code{\link{doubleFluxDel}} and \code{\link{doubleGeneDel}})
  and fills a solution object (a derivative of class
  \code{\linkS4class{optsol}}) and a problem object depending on the
  specified solver. The optimizations are then performed by calls to
  \code{\link{optimizeProb}}.
}

\usage{
optimizer(model, lb, ub,
          delete, geneFlag,
          algorithm = SYBIL_SETTINGS("ALGORITHM"),
          setToZero = FALSE,
          checkOptSolObj = FALSE,
          rebuildModel = FALSE,
          fld = "none",
          prCmd = NA, poCmd = NA,
          prDIR = NULL, poDIR = NULL,
          verboseMode = 2,
          ...)
}

\arguments{
  \item{model}{
    An object of class \code{modelorg}.
  }
  \item{lb}{
    A numeric vector containing the lower bounds for the reactions/genes in the
    \code{dels} slot of \code{\linkS4class{optsol}}.  If \code{lb} is
    missing, the values are set to 0.
  }
  \item{ub}{
    A numeric vector containing the upper bounds for the reactions/genes in the
    \code{dels} slot of \code{\linkS4class{optsol}}.  If \code{ub} is
    missing, the values are set to 0.
  }
  \item{delete}{
    A matrix containing gene id's or indices of fluxes to delete. The deletions
    were performed line by line.
  }
  \item{geneFlag}{
    A boolean flag. If set to \code{TRUE}, all entries in \code{delete} are
    considered to be gene id's, otherwise, they were treated as fluxes.
  }
  \item{fld}{
    Type of flux distribution to return. If set to \code{"none"}, no flux
    distribution will be returned. If set to \code{"fluxes"}, only the real
    flux distribution is returned, meaning all variable values after
    optimization representing a flux (reaction) in the model. If set to
    \code{"all"}, all variable values are returned.
  }
  \item{algorithm}{
    A single character value giving the algorithm to compute genetic
    pertubations. Can be
    \code{"\link[=sysBiolAlg_fba-class]{fba}"}: flux-balance analysis,
    \code{"\link[=sysBiolAlg_moma-class]{moma}"}: minimization of
                                                  metabolic adjustment (MOMA),
    \code{"\link[=sysBiolAlg_lmoma-class]{lmoma}"}: linear version of MOMA or
    \code{"\link[=sysBiolAlg_room-class]{room}"}: regulatory on/off
                                                  minimization (ROOM).\cr
    Default: \code{SYBIL_SETTINGS("ALGORITHM")}.
  }
  \item{setToZero}{
    Logical: If the lp-solver returns a solution which is not optimal, set the
    corresponding objective value to zero.\cr
    Default: \code{FALSE}.
  }
  \item{verboseMode}{
    Single integer value, giving the amount of output to the console. Use
    \code{\link{sink}} to redirect output to a file. If \code{loglevel == 2},
    a progress bar will be produced. If \code{loglevel > 2}, intermediate
    results will be printed.\cr
    Default: \code{2}.
  }
  \item{checkOptSolObj}{
    Boolean: check the optsol object with \code{\link{checkOptSol}}.\cr
    Default: \code{FALSE}.
  }
  \item{rebuildModel}{
    Logical.  If set to \code{TRUE}, the problem object will be rebuilt prior
    each round of optimization.
    Default: \code{FALSE}.
  }
  \item{prCmd}{
    A list of preprocessing commands passed to \code{\link{optimizeProb}}.
    See there for details.\cr
    Default: \code{NA}.
  }
  \item{poCmd}{
    A list of postprocessing commands passed to \code{\link{optimizeProb}}.
    See there for details.\cr
    Default: \code{NA}.
  }
  \item{prDIR}{
    A numeric or character vector, indicating in which round of optimization
    the preprocessing command(s) will be executed. \code{prDIR = c(2, 5, 10)}
    executes the commands in prCmd before the second, 5th and 10th optimization.
    
    If \code{prDIR} is a character vector, for example
    \code{prDIR = c("10")}, the preprocessing commands given in \code{prCmd}
    will pe executed every 10th round of optimization. 
    
    If \code{prDIR} is character and has length \code{2}, the first element is
    an offset to the following elements. \code{prDIR = c("-2", "10")} will do
    the preprocessing on every 10th round of optimization, beginning in round
    number \code{10 - 2 = 8}.\cr
    Default: \code{NULL}.
  }
  \item{poDIR}{
    The same as \code{prDIR}, but for postprocessing.
    Default: \code{NULL}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{sysBiolAlg}}.
  }
}

\details{
  The function \code{optimizer} is called by the functions
  \code{\link{oneFluxDel}}, \code{\link{oneGeneDel}},
  \code{\link{doubleFluxDel}}, \code{\link{doubleGeneDel}} and
  \code{\link{geneDeletion}}, but not directly by a user. The function
  \code{optimizer} runs a series of optimizations, whereas the function
  \code{\link{optimizeProb}} only runs exactly one optimization. In principle,
  \code{optimizer} is a kind of wrapper to \code{\link{optimizeProb}}.
}

\value{
  A derivative of class \code{\linkS4class{optsol}}.
}

\references{
  Becker, S. A., Feist, A. M., Mo, M. L., Hannum, G., Palsson, B. Ø.  and
  Herrgard, M. J. (2007) Quantitative prediction of cellular metabolism with
  constraint-based models: the COBRA Toolbox. \emph{Nat Protoc} \bold{2},
  727--738.

  Schellenberger, J., Que, R., Fleming, R. M. T., Thiele, I., Orth, J. D.,
  Feist, A. M., Zielinski, D. C., Bordbar, A., Lewis, N. E., Rahmanian, S.,
  Kang, J., Hyduke, D. R. and Palsson, B. Ø. (2011) Quantitative prediction of
  cellular metabolism with constraint-based models: the COBRA Toolbox v2.0.
  \emph{Nat Protoc} \bold{6}, 1290--1307.
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\link{oneFluxDel}}, \code{\link{oneGeneDel}},
         \code{\link{doubleFluxDel}}, \code{\link{doubleGeneDel}},
         \code{\linkS4class{optsol}}, \code{\link{optimizeProb}}
         and \code{\link{SYBIL_SETTINGS}}.}


\keyword{optimize}

