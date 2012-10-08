\name{modelorg2ExPA}
\alias{modelorg2ExPA}

\encoding{utf8}

\title{
  Write an Instance of Class modelorg to File in ExPA Format
}

\description{
  The function \code{modelorg2ExPA} writes the content of an instance of class
  \code{\link{modelorg}} to text files in a format which can be read by the
  program ExPA to compute extreme pathways.
}

\usage{
  modelorg2ExPA(model, fname, exIntReact = NA, filepath = ".", suffix = "expa")
}

\arguments{
  \item{model}{
    An object of class \code{\link{modelorg}}.
  }
  \item{fname}{
    An single character string giving the filename to write to.\cr
    Default: \code{<model_id>.expa}.
  }
  \item{exIntReact}{
    An object of class \code{\link{reactId}}, character or integer, giving id's
    of internal reactions to exclude in the ExPA file.\cr
    Default: \code{NA}.
  }
  \item{filepath}{
    A single character string giving the path to a certain directory in which
    the output files will be stored.\cr
    Default: \code{"."}.
  }
  \item{suffix}{
    A single character string giving the file name extension.\cr
    Default: \code{"expa"}.
  }
}

\details{
  The function \code{modelorg2ExPA} produces input files for the program ExPA.
  With ExPA, it is possible to calculate extreme pathways in metabolic networks.
}

\value{
  Returns \code{TRUE} on success.
}

\references{
  Bell, S. L. and Palsson, B. Ø. (2005) Expa: a program for calculating extreme
  pathways in biochemical reaction networks. \emph{Bioinformatics} \bold{21},
  1739--1740.

  The ExPA homepage \url{http://gcrg.ucsd.edu/Downloads/ExtremePathwayAnalysis}.
}

\author{
  Gabriel Gelius-Dietrich, C. Jonathan Fritzemeier
}

\keyword{ IO }
