\name{geneDeletion}
\alias{geneDeletion}

\encoding{utf8}

\title{
  Simulate N Gene Deletion Experiments
}

\description{
  The function \code{geneDeletion} studies the effect of in sillico gene
  deletions on the current objective function (e.g. biomass production).
}

\usage{
  geneDeletion(model, genes, combinations = 1, ...)
}

\arguments{
  \item{model}{
    An object of class \code{\linkS4class{modelorg}}.
  }
  \item{genes}{
    Character or Integer: the genes to delete (see Details below).
  }
  \item{combinations}{
    A single integer value. If \code{combinations > 1} and \code{genes} is not a
    matrix, \code{combinations} is the number of elements from \code{genes}
    taken at a time while building all combinations of the elements in
    \code{genes} (see Details below).\cr
    Default: \code{1}.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{optimizer}}.
  }
}

\details{
  If argument \code{genes} is a matrix of character values (gene id's) or
  integers (pointers to gene id's), each column is treated as one deletion
  experiment. If the matrix is made up of integers, a zero entry means no gene.
  
  If argument \code{genes} is a character vector or integer, the argument
  \code{combinations} gives the number of gene id's taken each time in order
  to build all possible combinations of \code{genes}. A matrix is constructed
  using \code{\link{combn}}.
  
  If argument \code{genes} is empty, the number of unique genes in \code{model}
  is used.
}

\value{
  An object of class \code{\link[=optsol_genedel-class]{optsol_genedel}}.
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\linkS4class{modelorg}}, \code{\linkS4class{optsol}},
  \code{\link[=optsol_genedel-class]{optsol_genedel}},
  \code{\link{checkOptSol}}, \code{\link{oneGeneDel}},
  \code{\link{optimizer}}, \code{\link{optimizeProb}},
  \code{\link{combn}} and \code{\link{SYBIL_SETTINGS}}.
}

\examples{
  ## load the dataset
  data(Ec_core)
  
  ## perform a single gene deletion analysis (delete every gene one time)
  gd <- geneDeletion(Ec_core)
  
  ## triple gene deletion analysis using the first ten genes
  gd <- geneDeletion(Ec_core, genes = 10, combinations = 3)

\dontrun{
  ## perform a double gene deletion analysis (delete all possible pairwise
  ## combinations of all genes)
  gd <- geneDeletion(Ec_core, combinations = 2)

  ## perform a triple gene deletion analysis
  ## (very high number of optimizations)
  gd <- geneDeletion(Ec_core, combinations = 3)
}
}

\keyword{optimize}

