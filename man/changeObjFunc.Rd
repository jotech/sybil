\name{changeObjFunc}

\alias{changeObjFunc}

\encoding{utf8}

\title{Sets/changes the Objective Function}

\description{
  The function \code{changeObjFunc} changes or sets the objective
  function for a specified model.
}

\usage{
  changeObjFunc(lpmodel, react, obj_coef, checkIds = TRUE)
}

\arguments{
  \item{lpmodel}{A model. An object of class
    \code{\linkS4class{modelorg}}, or a problem object of a lp solver.}
  \item{react}{Charater vector containing reaction id's, or a numerical
    vector containing indices of reaction id's.}
  \item{obj_coef}{A numerical vector with the same length as
    \code{react} containing the objective coefficients. If missing, a
    vector of the same length as \code{react} with all ones (1) is produced.}
  \item{checkIds}{Boolean: check \code{react} if ok (default: TRUE)}
}

\details{
  The argument \code{react} will be evaluated by the function
  \code{\link{checkReactId}}. The return value is used to change
  the objective function.
}

\value{
  Returns the given model (an object of the same class as the argument
  \code{lpmodel}) containing the new objective function.
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\link{checkReactId}}}

\examples{
  ## sets the objective function to the ATP maintenance reaction:
  data(Ec_core)
  Ec_new <- changeObjFunc(Ec_core, "ATPM")
}

\keyword{optimize}

