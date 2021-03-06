\name{optObj_glpkAPI-class}

\Rdversion{1.1}
\encoding{utf8}
\docType{class}

\alias{optObj_glpkAPI-class}

\title{Class \code{"optObj_glpkAPI"}}

\description{
  Structure of the class \code{"optObj_glpkAPI"}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the constructor function \code{optObj}:

  \code{test <- optObj(solver = "glpkAPI")}.
}

\section{Slots}{
  \describe{
    \item{\code{oobj}:}{
      Object of class \code{"pointerToProb"} containing a pointer to a
      \pkg{glpkAPI} problem object.
    }
    \item{\code{solver}:}{
      Object of class \code{"character"} containing the name of the solver
      software (see \code{\link{SYBIL_SETTINGS}} for suitable values).
    }
    \item{\code{method}:}{
      Object of class \code{"character"} containing the method (algorithm) used
      by the solver software (see \code{\link{SYBIL_SETTINGS}} for suitable
      values).
    }
    \item{\code{probType}:}{
      Object of class \code{"character"} giving the problem type (see
      \code{\link{optObj}} for suitable values).
    }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{optObj}"}, directly.
}

\author{
  Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>

  Maintainer: Claus Jonathan Fritzemeier <clausjonathan.fritzemeier@uni-duesseldorf.de>
}


\seealso{
  Superclass \code{\linkS4class{optObj}} and constructor function
  \code{\link{optObj}}
}

\examples{
  showClass("optObj_glpkAPI")
}

\keyword{classes}


