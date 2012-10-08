\name{geneDel}
\alias{geneDel}

\encoding{utf8}

\title{
  Get Gene-Reaction Association
}

\description{
  The function \code{geneDel} finds the fluxes corresponding to gene id's.
}

\usage{
  geneDel(model, genes, checkId = FALSE)
}

\arguments{
  \item{model}{
    An object of class \code{modelorg}.
  }
  \item{genes}{
    A vector of character strings of gene id's used in \code{model}, or an
    integer vector with indices to gene id's in \code{\link{allGenes}(model)}.
  }
  \item{checkId}{
    Boolean. If set to \code{TRUE}, argument \code{genes} will be checked wether
    it fits to \code{model} (e.g. are all genes existing). If set to
    \code{FALSE}, \code{genes} must contain indices of gene id's in
    \code{model}, e.g. in calls from \code{\link{optimizer}}.
  }
}

\details{
  The function \code{geneDel} checks for a set of gene id's in \code{gene}
  which fluxes this set effects.
}

\value{
  An integer vector of pointers to reaction id's in \code{model} or \code{NA}.
}

\references{
  Edwards, J. S., Ibarra, R. U. and Palsson, B. Ø. (2001) In silico predictions
  of \emph{Escherichia coli} metabolic capabilities are consistent with
  experimental data. \emph{Nat Biotechnol} \bold{19}, 125--130.
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  \code{\link{optimizer}}
}

