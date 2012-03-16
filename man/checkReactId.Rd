\name{checkReactId}

\alias{checkReactId}

\encoding{utf8}

\title{Check if a Reaction Id is Valid}

\description{
  The function \code{checkReactId} evaluates a vector of reaction id's
  if they are unique and appear in a given model.
}

\usage{
  checkReactId(model, react, needId = FALSE)
}

\arguments{
  \item{model}{A model. An object of class
    \code{\linkS4class{modelorg}}, or a problem object of a lp solver.}
  \item{react}{Charater vector containing reaction id's, or a numerical
    vector containing indices of reaction id's.}
  \item{needId}{Logical: if reaction ids are needed (default: FALSE).}
}

\details{
  If argument \code{react} is numeric, the maximun value will be
  inspected, if it is larger than the number of reactions in the model.

  In case of a character vector, \code{react} is matched to the reaction
  id's residing in the model. If they are not found, grep is used.

  If argument \code{react} is of class \code{\linkS4class{reactId}}, it will
  be returned without checking.
}

\value{
  Returns NA if the values of \code{react} do not correspond to the
  model, or an object of class \code{\linkS4class{reactId}}.
}

\author{Gabriel Gelius-Dietrich}

\note{
  If the found reaction id's are not unique, NA is returned.
}

\seealso{\code{\linkS4class{reactId}}}

\examples{
  data(Ec_core)
 
  ## Example with react as character vector
  ids <- c("ATPM", "ACK")
  idc <- checkReactId(Ec_core, ids)

  ## Example with react as numerical vector
  ids <- c(1:4)
  idc <- checkReactId(Ec_core, ids)
}

\keyword{character}

