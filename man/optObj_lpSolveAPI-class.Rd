\name{optObj_lpSolveAPI-class}
\Rdversion{1.1}
\docType{class}
\alias{optObj_lpSolveAPI-class}

\alias{return_codeLPSOLVE}
\alias{loadMatrixPerColumnLPSOLVE}

\title{Class \code{"optObj_lpSolveAPI"}}

\description{
  Structure of the class \code{"optObj"}. Objects of that class are
  returned by the function \code{\link{prepProbObj}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{optObj}:

  \code{test <- optObj(solver = "lpSolveAPI")}.
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

\section{Extends}{
  Class \code{"\linkS4class{optObj}"}, directly.
}

\section{Methods}{
  No methods defined with class "optObj_lpSolveAPI" in the signature.
}

\section{Further usefull Functions}{
  \describe{
    \item{\code{return_codeLPSOLVE}:}{(\code{code})
        prints a human readable translation of return codes of lpSolveAPI. }
    \item{\code{loadMatrixPerColumnLPSOLVE}:}{(\code{lpmod, constMat})
        load a constraint matrix to lpSolveAPI column by column. }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\link{prepProbObj}}}

\examples{
  showClass("optObj_lpSolveAPI")
}

\keyword{classes}


