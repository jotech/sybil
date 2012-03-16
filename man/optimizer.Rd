\name{optimizer}
\alias{optimizer}

\encoding{utf8}

\title{Performs Optimizations for the Gene/Flux Deletion Experiments}

\description{
  This function is used by the gene/flux deletion experiments
  (the functions \code{\link{oneFluxDel}}, \code{\link{oneGeneDel}},
  \code{\link{doubleFluxDel}} and \code{\link{doubleGeneDel}})
  and fills a solution object (a derivative of class \code{\link{optsol}})
  and a problem object depending on the specified solver.
  The optimizations are then performed by calls to \code{\link{simpleFBA}}.
}

\usage{
optimizer(model, optsol, lb, ub,
          wtFluxes = NA,
          alg = SYBIL_SETTINGS("ALGORITHM"),
          setToZero = FALSE,
          checkOptSolObj = FALSE,
          rebuildModel = FALSE,
          copyModel = FALSE,
          prCmd = NA, poCmd = NA,
          prDIR = NA, poDIR = NA,
          verboseMode = 3,
          loglevel = -1,
          logfile = NA,
          logfileEnc = NA,
          ...)
}

\arguments{
  \item{model}{
    An object of class \code{modelorg}.
  }
  \item{optsol}{
    An object of a class which inherits class \code{\linkS4class{optsol}}.
  }
  \item{lb}{
    A numeric vector containing the lower bounds for the reactions/genes in the
    \code{dels} slot of \code{optsol}.  If \code{lb} is missing, the values are
    set to 0.
  }
  \item{ub}{
    A numeric vector containing the upper bounds for the reactions/genes in the
    \code{dels} slot of \code{optsol}.  If \code{ub} is missing, the values are
    set to 0.
  }
  \item{wtFluxes}{
    A single numeric value, or a numeric vector of length equal to the number of
    fluxes (variables) in the model (the argument has the same meaning
    as in \code{\link{simpleFBA}}).\cr
    Default: \code{NA}.
  }
  \item{alg}{
    A single character value. Can be \dQuote{\code{FBA}},
    \dQuote{\code{linearMOMA}} or \dQuote{\code{linearMOMA_COBRA}} (see Details
    below).\cr
    Default: \code{SYBIL_SETTINGS("ALGORITHM")}.
  }
  \item{setToZero}{
    Logical: If the lp-solver returns a solution which is not optimal, set the
    corresponding objective value to zero.\cr
    Default: \code{FALSE}.
  }
  \item{verboseMode}{
    Single integer value, see \code{\link{sybilLog}} for possible values.
    In addition, \code{loglevel > 2} prints a progress bar.\cr
    Default: \code{3}.
  }
  \item{loglevel}{
    Single integer value, see \code{\link{sybilLog}} for possible values.
    Default: \code{-1}.
  }
  \item{logfile}{
    A single character string giving the filename of the logfile.
    If \code{logfile} is set to \code{NA} and \code{loglevel} is \code{> -1},
    a default filename is used (see \code{\link{sybilLog}} for details).\cr
    Default: \code{NA}.
  }
  \item{logfileEnc}{
    A single character string giving the encoding for the logfile. If set to
    \code{NA}, the default encoding is used: \code{getOption("encoding")}.
    If \code{logfileEnc} is set to something different, the option
    \code{useFancyQuotes} is turned
    off: \code{options(useFancyQuotes = FALSE)}.\cr
    Default: \code{NA}.
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
  \item{copyModel}{
    Logical.  If set to \code{TRUE}, the problem object will be copied from a
    \dQuote{backup copy} of the initial problem object prior each round of
    optimization.
    Default: \code{FALSE}.
  }
  \item{prCmd}{
    A list of preprocessing commands passed to \code{\link{simpleFBA}}.
    See there for details.\cr
    Default: \code{NA}.
  }
  \item{poCmd}{
    A list of postprocessing commands passed to \code{\link{simpleFBA}}.
    See there for details.\cr
    Default: \code{NA}.
  }
  \item{prDIR}{
    A numeric or character vector, indicating in which round of optimization
    the preprocessing command(s) will be executed. \code{prDIR = c(2, 5, 10)}
    executes the commands in prCmd before the second, 5th and 10th optimization.
    
    If \code{prDIR} is a character vector, for example
    \code{prDIR = c("10")} executes the commands in \code{prCmd} every 10th
    round of optimization.  If the last element of \code{prDIR} is \code{"WT"},
    preprocessing is additionally performed on the wild type (first) problem.
    If \code{prDIR} is character and has a length \code{> 1} (a possible last
    element \code{"WT"} excluded), the first element is an offset to the
    following elements. \code{prDIR = c("-2", "10")} will do the preprocessing
    on every 10th round of optimization, beginning in round number
    \code{10 - 2 = 8}.\cr
    Default: \code{NA}.
  }
  \item{poDIR}{
    The same as \code{prDIR}, but for postprocessing.
    Default: \code{NA}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{prepProbObj}}.
  }
}

\details{
  The function \code{optimizer} is called by the functions
  \code{\link{oneFluxDel}}, \code{\link{oneGeneDel}},
  \code{\link{doubleFluxDel}} and \code{\link{doubleGeneDel}},
  but not by a user.  The important arguments here are \code{resolve},
  \code{setToZero}, \code{logfile} and \code{verboseMode} which can
  be set by using \code{\dots} in the invoking functions.
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
         \code{\linkS4class{optsol}}, \code{\link{prepProbObj}}
         and \code{\link{SYBIL_SETTINGS}}.}


\keyword{optimize}

