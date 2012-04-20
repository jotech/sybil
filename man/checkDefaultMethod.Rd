\name{checkDefaultMethod}
\alias{checkDefaultMethod}

\encoding{utf8}

\title{
  Validate Solver and Method
}

\description{
  The function \code{checkDefaultMethod} returns the default method
  for a desired solver, or a default solver -- method pair.
}

\usage{
  checkDefaultMethod(solver, method)
}

\arguments{
  \item{solver}{
    A single character value, containing the solver name
    (see \code{\link{SYBIL_SETTINGS}}).
  }
  \item{method}{
    A single character value, containing the method name
    (see \code{\link{SYBIL_SETTINGS}}).
  }
}

\details{
  At first, the function checks if argument \code{solver} contains a valid
  solver, if not the default solver will be returned. Next the same is done
  for the argument \code{method}, regarding the current value of \code{solver}.
}

\value{
  \item{sol}{Validated solver name.}
  \item{met}{Validated method name.}
  \item{parm}{Default prarmeter set for the validated method.}
}

\author{
  Gabriel Gelius-Dietrich

  Maintainer: Gabriel Gelius-Dietrich <geliudie@uni-duesseldorf.de>
}

\seealso{
  \code{\link{SYBIL_SETTINGS}}
}

\keyword{ optimize }
