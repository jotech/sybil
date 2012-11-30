\name{sysBiolAlg_moma-class}

\Rdversion{1.1}
\encoding{utf8}

\docType{class}

\alias{sysBiolAlg_moma-class}
\alias{sysBiolAlg_moma}
\alias{moma}

\title{Class \code{"sysBiolAlg_moma"}}

\description{
  The class \code{sysBiolAlg_moma} holds an object of class
  \code{\linkS4class{optObj}} which is generated to meet the
  requirements of the MOMA algorithm.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  
  \code{sysBiolAlg(model, algorithm = "moma", ...)}.
  
  Arguments to \code{...} which are passed to method \code{initialize} of class
  \code{sysBiolAlg_moma} are described in the Details section.
}

\section{Slots}{
  \describe{
    \item{\code{problem}:}{Object of class \code{"optObj"}
      containing the problem object.
    }
    \item{\code{algorithm}:}{Object of class \code{"character"}
      giving the name of the algorithm.
    }
    \item{\code{nr}:}{Object of class \code{"integer"}
      giving the number of rows of the problem object.
    }
    \item{\code{nc}:}{Object of class \code{"integer"}
      giving the number of columns of the problem object
    }
    \item{\code{fldind}:}{Object of class \code{"integer"}
      pointers to columns (variables) representing a flux (reaction) in the
      original network. The variable \code{fldind[i]} in the problem object
      represents reaction \code{i} in the original network. 
    }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{sysBiolAlg}"}, directly.
}

\section{Methods}{
  No methods defined with class "sysBiolAlg_moma" in the signature.
}

\details{
  The \code{initialize} method has the following arguments:
  \describe{
    \item{model}{
      An object of class \code{\linkS4class{modelorg}}.
    }
    \item{wtflux}{
      A numeric vector holding an optimal wild type flux distribution for the
      given model. If missing, a default value is computed based on FBA.
      If given, arguments \code{solver} and \code{method} are used, but
      \code{solverParm} is not.
    }
    \item{scaling}{
      Scaling options used to scale the constraint matrix. If set to
      \code{NULL}, no scaling will be performed
      (see \code{\link{scaleProb}}).\cr
      Default: \code{NULL}.
    }
    \item{...}{
      Further arguments passed to the initialize method of
      \code{\linkS4class{sysBiolAlg}}. They are \code{solver},
      \code{method} and \code{solverParm}.
    }
  }

  The problem object is built to be capable to perform the MOMA algorithm with
  a given model, which is basically the solution of a quadratic programming
  problem
  \deqn{%
    \begin{array}{rll}%
      \min            & \begin{minipage}[b]{5em}
                        \[ 
                          \sum_{i,j=1}^n
                          \bigl(v_{j,\mathrm{del}} - v_{i,\mathrm{wt}}\bigr)^2
                        \]
                        \end{minipage}                                 \\[2em]
          \mathrm{s.\,t.} & \mbox{\boldmath$Sv$\unboldmath} = 0         \\[1ex]
                          & \alpha_i \leq v_i \leq \beta_i
                            & \quad \forall i \in \{1, \ldots, n\}      \\[1ex]
    \end{array}%
  }{
     min  sum (v_j,del - v_i,wt)^2  for i,j = 1, ..., n
      s.t. Sv = 0
           a_i <= v_i <= b_i  for i = 1, ..., n
  }
  with \eqn{\bold{S}}{S} beeing the stoichiometric matrix, \eqn{\alpha_i}{a_i}
  and \eqn{\beta_i}{b_i} beeing the lower and upper bounds for flux (variable)
  \eqn{i} (\eqn{j} for the deletion strain).
  The total number of variables of the optimization problem is denoted
  by \eqn{n}.
  Here,
  \eqn{
    \mbox{\boldmath$v$\unboldmath}_{\mathrm{wt}}
  }{v_wt}
  is the optimal wild type flux distribution. This can be set via the argument
  \code{wtflux}.  If \code{wtflux} is \code{NULL} (the default), the
  wild type flux distribution will be calculated by a standard FBA.
  The optimization can be executed by using \code{\link{optimizeProb}}.
}


\references{
  Segrè, D., Vitkup, D. and Church, G. M. (2002) Analysis or optimality in
  natural and pertubed metabolic networks. \emph{PNAS} \bold{99}, 15112--15117.
}

\author{
  Gabriel Gelius-Dietrich
}

\seealso{
  Constructor function \code{\link{sysBiolAlg}} and
  superclass \code{\linkS4class{sysBiolAlg}}.
}

\examples{
  showClass("sysBiolAlg_moma")
}

\keyword{classes}
