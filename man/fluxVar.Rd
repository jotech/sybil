\name{fluxVar}
\alias{fluxVar}

\encoding{utf8}

\title{
  Flux Variability Analysis
}
\description{
  Performs flux variablilty analysis for a given model.
}
\usage{
  fluxVar(model, react, fld = FALSE, verboseMode = 2, ...)
}

\arguments{
  \item{model}{
    An object of class \code{\link{modelorg}}.
  }
  \item{react}{
    An object of class \code{\linkS4class{reactId}}, character or integer.
    Specifies the fluxes (variables) to analyse.\cr
    If left empty: \code{react_id(model)}.
  }
  \item{fld}{
    Boolean.  Save the resulting flux distribution.\cr
    Default: \code{FALSE}
  }
  \item{verboseMode}{
    An integer value indicating the amount of output to stdout:
    0: nothing, 1: status messages, 2: like 1 plus a progress indicator,
    3: a table containing the reaction id's and the corresponding
    min max values.\cr
    Default: \code{2}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{sysBiolAlg}}.
    See class \code{\link[=sysBiolAlg_fv-class]{sysBiolAlg_fv}} for details.
  }
}

\details{
  The algorithm is described in
  \code{\link[=sysBiolAlg_fv-class]{sysBiolAlg_fv}}.
}

\value{
  An object of class \code{\link[=optsol_fluxVar-class]{optsol_fluxVar}}.
}

\references{
  Becker, S. A., Feist, A. M., Mo, M. L., Hannum, G., Palsson, B. Ø.  and
  Herrgard, M. J. (2007) Quantitative prediction of cellular metabolism with
  constraint-based models: the COBRA Toolbox. \emph{Nat Protoc} \bold{2},
  727--738.

  Schellenberger, J., Que, R., Fleming, R. M. T., Thiele, I., Orth, J. D.,
  Feist, A. M., Zielinski, D. C., Bordbar, A., Lewis, N. E., Rahmanian, S.,
  Kang, J., Hyduke, D. R. and Palsson, B. Ø. (2011) Quantitative prediction of
  cellular metabolism with constraint-based models: the COBRA Toolbox v2.0.
  \emph{Nat Protoc} \bold{6}, 1290--1307.

  Bernhard Ø. Palsson (2006).
  \emph{Systems Biology: Properties of Reconstructed Networks}.
  Cambridge University Press.
}

\author{Gabriel Gelius-Dietrich}

\examples{
  data(Ec_core)
  fv <- fluxVar(Ec_core)
  plot(fv)
}

\keyword{optimize}

