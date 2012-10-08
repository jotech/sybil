\name{checkOptSol}

\alias{checkOptSol}

\encoding{utf8}

\title{Interprets an Object of Class Optsol}

\description{
  The function \code{checkOptSol} evaluates the results of the solution
  of a linear optimization; the returned objects e.g. from
  \code{\link{optimizeProb}}.
}

\usage{
  checkOptSol(optsol, onlywarn = FALSE)
}

\arguments{
  \item{optsol}{An object of class \code{\link[=optsol-class]{optsol}}.}
  \item{onlywarn}{Logical, default: FALSE.}
}

\details{
  The function \code{checkOptSol} is used by functions performing a
  linear optimization (e.g. \code{\link{optimizeProb}}). In that case, the
  argument \code{onlywarn} is set to TRUE. If the optimization ends
  unsuccesfull, a warning will be produced.

  It is also possible to use the function directly, with \code{onlywarn}
  set to FALSE (the default). In that case, an object of class
  \code{\linkS4class{checksol}} will be retuned. This object contains a summary
  with the exit status of the optimization.
}

\value{
  TRUE or FALSE if \code{onlywarn} is set to TRUE, an object of class
  \code{\linkS4class{checksol}} otherwise.
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\linkS4class{checksol}},
  \code{\link{optimizeProb}} and \code{\link{oneFluxDel}}
}

\examples{
  data(Ec_core)
  Ec_f <- optimizeProb(Ec_core, retOptSol = TRUE)
  Ec_check <- checkOptSol(Ec_f)
}

\keyword{optimize}

