\name{sybil-package}

\alias{sybil-package}
\alias{sybil}

\encoding{utf8}

\docType{package}

\title{
  Systems Biology Library
}

\description{
  The package \code{sybil} is a collection of functions designed for
  in silico analysis of metabolic networks.
}

\details{

  The package SyBiL is designed to read csv files. This is done by
  the function \code{readTSVmod}. This function returns
  an object of the class \code{modelorg}.

  the beginning (or reading a csv file):\cr
  \code{ec <- readTSVmod(reactList = "Ec_core.csv")}

  find double reactions:\cr
  \code{ec_d <- doubleReact(ec)}

  remove reactions one and five:\cr
  \code{ec_r <- rmReact(ec, c(1, 5))}

  convert to irreversible model:\cr
  \code{ec_ir <- mod2irrev(ec)}

  perform FBA:\cr
  \code{ec_f <- simpleFBA(ec)}

  perform single gene deletion analysis:\cr
  \code{ec_o <- oneGeneDel(ec)}

}

\author{
  Gabriel Gelius-Dietrich

  Maintainer: geliudie@uni-duesseldorf.de
}

\references{
  The BiGG database \url{http://bigg.ucsd.edu/}.

  Schellenberger, J., Park, J. O., Conrad, T. C., and Palsson, B. Ø., (2010)
  BiGG: a Biochemical Genetic and Genomic knowledgebase of large scale metabolic
  reconstructions. \emph{BMC Bioinformatics} \bold{11}, 213.

  Becker, S. A., Feist, A. M., Mo, M. L., Hannum, G., Palsson, B. Ø.  and
  Herrgard, M. J. (2007) Quantitative prediction of cellular metabolism with
  constraint-based models: the COBRA Toolbox. \emph{Nat Protoc} \bold{2},
  727--738.

  Schellenberger, J., Que, R., Fleming, R. M. T., Thiele, I., Orth, J. D.,
  Feist, A. M., Zielinski, D. C., Bordbar, A., Lewis, N. E., Rahmanian, S.,
  Kang, J., Hyduke, D. R. and Palsson, B. Ø. (2011) Quantitative prediction of
  cellular metabolism with constraint-based models: the COBRA Toolbox v2.0.
  \emph{Nat Protoc} \bold{6}, 1290--1307.
}

\seealso{
  Package \code{sybilSBML} (\code{readSBMLmod}).% (\code{\link[sybilSBML]{readSBMLmod}}).
}

\examples{
  data(Ec_core)
  Ec_ofd <- oneFluxDel(Ec_core)
  hist(lp_obj(Ec_ofd))
}

\keyword{package}

