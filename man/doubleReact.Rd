\name{doubleReact}
\alias{doubleReact}

\encoding{utf8}

\title{Identifies Identical Reactions}

\description{
  The function \code{doubleReact} identifies identical reactions
  (isoenzymes) in a model.
}

\usage{
  doubleReact(model, checkRev = TRUE)
}

\arguments{
  \item{model}{An object of class \code{\linkS4class{modelorg}}.}
  \item{checkRev}{Logical: Check reversibility of reactions (default: TRUE).}
}

\details{
  In the first step, the stoichiometric matrix S is devided into
  groups of reactions containing the same number of metabolites. After
  that, the row indices of the non-zero elements of these matrixces are
  compared. If identical pairs are found, we check the corresponding
  values in S. If they are also identical, the reversibility of the reactions
  are examined. If they are the same, the two reactions are called identical.
}

\value{
  If no identical reactions were found, the return value is
  FALSE. Otherwise a list is returned, ordered by the number of
  metabolites used in each reaction. The first element is the number of
  metabolites. The following elements are one or more numerical vectors
  containing the indices of the reactions.
}


\author{Gabriel Gelius-Dietrich}

\note{
  At the moment, the directions of a pair of reactions is not
  compared. Meaning, that if concerning to the values in S the reaction
  is in forward direction, but not when including the flux values,
  \code{doubleReact} will not find it.

  There are some things left to do with the upper and lower bounds.
}

\examples{
  data(Ec_core)
  Ec_d <- doubleReact(Ec_core)
}

\keyword{manip}
