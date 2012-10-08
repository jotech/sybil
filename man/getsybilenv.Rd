\name{getsybilenv}
\alias{getsybilenv}

\title{
  Print sybil Environment
}

\description{
  Prints current settings in the sybil environment.
}

\usage{
  getsybilenv(part)
}

\arguments{
  \item{part}{
    A character vector containing names of elements in the sybil environment.
    Possible values are:
    \describe{
      \item{\code{"solvers"}}{
        supported R packages for solving optimization problems.
      }
      \item{\code{"methods"}}{
        methods to solve optimization problems included in the R packages.
      }
      \item{\code{"ptype"}}{
        methods required for a particular problem type.
      }
      \item{\code{"algorithm"}}{
        algorithms used in systems biology to use with particular functions.
      }
    }
  }
}

\details{
  Typical usages are
  \preformatted{
    getsybilenv(PART)
    getsybilenv()
  }
  If argument \code{part} is not given, all four elements described above where
  printed.
}

\value{
  Returns \code{NULL} invisibly.
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{addAlgorithm}}, \code{\link{addSolver}},
  \code{\link{SYBIL_SETTINGS}}, \code{\link{checkDefaultMethod}} and
  \code{\link{checkAlgorithm}}
}

