\name{doubleGeneDel}
\alias{doubleGeneDel}

\encoding{utf8}

\title{Double Gene Deletion Experiment}

\description{
  Sets every pairwise combination of genes in a model for one time to zero.
}

\usage{
  doubleGeneDel(model, geneList1, geneList2, allComb = FALSE,
                lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                solver = SYBIL_SETTINGS("SOLVER"),
                method = SYBIL_SETTINGS("METHOD"),
                fld = FALSE, exLethal = TRUE,
                tol = SYBIL_SETTINGS("TOLERANCE"), ...)
}

\arguments{
  \item{model}{An object of class \code{modelorg}.}
  \item{geneList1}{Character: the genes to delete (default: \code{allGenes(model)}).}
  \item{geneList2}{Character: the genes to delete (default: \code{allGenes(model)}).}
  \item{allComb}{Logical: If set to \code{TRUE}, all possible pairwise combinations of geneList1 and geneList2 were deleted.
    If set to \code{FALSE}, geneList1 and geneList2 must have the same length. Every pair c(geneList1[i], geneList2[i]) will be deleted (default: FALSE).}
  \item{lpdir}{Direction of optimization (default: max).}
  \item{solver}{The linear problem solver (default: glpk).}
  \item{method}{The method for the linear problem solver (default: simplex).}
  \item{fld}{Logical: If set to TRUE, the flux distribution for
    each optimization step will be stored. This increases the running
    time significantly (default: FALSE).}
  \item{exLethal}{Logical: If set to TRUE, lethal genes are removed from the
    analysis. This increases the running time significantly (default: TRUE).}
  \item{tol}{Numeric: threshold value for a gene beeing lethal or not (default: 1E-6).}
  \item{...}{Further arguments passed to \code{\link{optimizer}}.}
}

\details{
  The function \code{doubleGeneDel} studies the effect of deletion of pairs of genes on
  the current objective function (e.g. biomass production).
  
  The first step is the calculation of a reference solution, which is
  the wild type solution. Then n further optimisations are performed (n:
  number of required optimizations) with \eqn{v_i = 0} for \eqn{i = 1,\dots,n}.
}

\value{
  An object of class \code{\link{optsol_doublegenedel}}.
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\link{modelorg}}, \code{\linkS4class{optsol}}, \code{\linkS4class{optsol_doublegenedel}}, \code{\link{checkOptSol}}, \code{\link{optimizer}} and  \code{\link{SYBIL_SETTINGS}}. }

\examples{
  data(Ec_core)
  Ec_ofd <- doubleGeneDel(Ec_core)
}

\keyword{optimize}

