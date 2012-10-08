\name{doubleGeneDel}
\alias{doubleGeneDel}

\encoding{utf8}

\title{Double Gene Deletion Experiment}

\description{
  Compute in sillico simulations of double-gene knock out mutants.
}

\usage{
  doubleGeneDel(model, geneList1, geneList2, allComb = FALSE,
                exLethal = TRUE, tol = SYBIL_SETTINGS("TOLERANCE"), ...)
}

\arguments{
  \item{model}{An object of class \code{modelorg}.}
  \item{geneList1}{Character: the genes to delete (default: \code{allGenes(model)}).}
  \item{geneList2}{Character: the genes to delete (default: \code{allGenes(model)}).}
  \item{allComb}{Logical: If set to \code{TRUE}, all possible pairwise combinations of geneList1 and geneList2 were deleted.
    If set to \code{FALSE}, geneList1 and geneList2 must have the same length. Every pair c(geneList1[i], geneList2[i]) will be deleted (default: FALSE).}
  \item{exLethal}{Logical: If set to TRUE, lethal genes are removed from the
    analysis. This increases the running time significantly (default: TRUE).}
  \item{tol}{Numeric: threshold value for a gene beeing lethal or not (default: 1E-6).}
  \item{\dots}{Further arguments passed to \code{\link{optimizer}}.}
}

\details{
  The function \code{doubleGeneDel} studies the effect of deletion of pairs of genes on
  the current objective function (e.g. biomass production).
  
  The first step is the calculation of a reference solution, which is
  the wild type solution. Then n further optimizations are performed (n:
  number of required optimizations) with \eqn{v_i = 0} for \eqn{i = 1,\dots,n}.
}

\value{
  An object of class \code{\link[=optsol_genedel-class]{optsol_genedel}}.
}

\author{Gabriel Gelius-Dietrich}

\seealso{
  \code{\link{modelorg}}, \code{\linkS4class{optsol}},
  \code{\link[=optsol_genedel-class]{optsol_genedel}},
  \code{\link{checkOptSol}}, \code{\link{optimizer}} and
  \code{\link{SYBIL_SETTINGS}}.
}

\examples{
  data(Ec_core)
  Ec_ofd <- doubleGeneDel(Ec_core)
}

\keyword{optimize}

