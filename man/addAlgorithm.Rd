\name{addAlgorithm}
\alias{addAlgorithm}

\title{
  Add a New Algorithm Name to sybil
}
\description{
  Certain simulations can be run using different algorithms. For example,
  a single gene deletion study (\code{\link{oneGeneDel}}) can be performed with
  FBA, MOMA or the like. With this funktion you can add a new algorithm to
  an existing kind of simulation.
}
\usage{
  addAlgorithm(alg, fkt = "optimizer")
}

\arguments{
  \item{alg}{
    A single character string containing the name of the new algorithm.
  }
  \item{fkt}{
    A function name, which can make use of the new algorithm.
  }
}

\value{
  Returns \code{NULL} invisibly.
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  \code{\link{checkAlgorithm}}, \code{\link{getsybilenv}}
}

