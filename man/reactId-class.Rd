\name{reactId-class}
\docType{class}

\alias{reactId-class}
\alias{reactId}
\alias{react_pos<-,reactId-method}
\alias{react_pos,reactId-method}
\alias{react_pos<-}
\alias{react_pos}
\alias{react_id<-,reactId-method}
\alias{react_id,reactId-method}
\alias{length,reactId-method}
\alias{show,reactId-method}

\encoding{utf8}

\title{Strukture of Class "reactId"}

\description{
  Structure of the class \code{"reactId"}. Objects of that class are
  returned by the function \code{\link{checkReactId}}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the function \code{reactId}:

  \code{test <- reactId(pos = c(1, 2), react = c("foo", "bar"))}.
}

\section{Slots}{
	 \describe{
    \item{\code{react_pos}:}{Object of class \code{"integer"} containing
      the index of reaction id's.}
    \item{\code{react_id}:}{Object of class \code{"character"}
      containing the reaction id's.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{react_pos<-}:}{\code{signature(object = "reactId")}:
      sets the \code{react_pos} slot.
    }
    \item{\code{react_pos}:}{
      \code{signature(object = "reactId")}: gets the \code{react_pos} slot.
    }
    \item{\code{react_id<-}:}{
      \code{signature(object = "reactId")}: sets the \code{react_id} slot.
    }
    \item{\code{react_id}:}{
      \code{signature(object = "reactId")}: gets the \code{react_id} slot.
    }
    \item{\code{length}}{
      \code{signature(object = "reactId")}: returns the number of reaction id's.
    }
    \item{\code{show}}{
      \code{signature(object = "reactId")}: prints the instance of class
      \code{reactId} as table.
    }
  }
}

\author{Gabriel Gelius-Dietrich}

\seealso{\code{\link{checkReactId}}}

\examples{
  showClass("reactId")
}

\keyword{classes}

