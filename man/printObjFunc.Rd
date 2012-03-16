\name{printObjFunc}
\alias{printObjFunc}

\encoding{utf8}

\title{Prints a Human Readable Objective Function}

\description{
  The function \code{printObjFunc} prints the objective function of a
  model in a human readable format.
}

\usage{
  printObjFunc(model)
}

\arguments{
  \item{model}{An object of class \code{\linkS4class{modelorg}}}
}

%\details{
%
%}

\value{
  Returns a character vector containing the objective function.
}


\author{Gabriel Gelius-Dietrich}

%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\seealso{\code{\link{changeObjFunc}}}

\examples{
  data(Ec_core)
  printObjFunc(Ec_core)

  ## with a different objective function than the default one
  ## just to illustrate how the function works ;-)
  Ec_objf <- changeObjFunc(Ec_core, c("EX_etoh(e)", "ETOHt2r"), c(1, 2))
  printObjFunc(Ec_objf)
}

\keyword{print}
