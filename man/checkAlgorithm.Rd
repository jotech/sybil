\name{checkAlgorithm}
\alias{checkAlgorithm}

\title{
  Check Algorithm
}

\description{
  Test, if a given algorithm can be used with a certain algorithm.
}

\usage{
  checkAlgorithm(alg, fkt = "optimizer")
}

\arguments{
  \item{alg}{
    A single character string containing the name of the algorithm.
  }
  \item{fkt}{
    A function name to test, if it can make use of the algorithm.
  }
}

\value{
  Returns \code{TRUE} if successfull, otherwise \code{FALSE}.
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  \code{\link{addAlgorithm}}, \code{\link{getsybilenv}}
}

