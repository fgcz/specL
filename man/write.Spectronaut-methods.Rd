\name{write.spectronaut-methods}

\docType{methods}

\alias{write.spectronaut}
\alias{write.spectronaut,specL-methods}
\alias{write.spectronaut,specLSet-methods}
\alias{write.spectronaut,ANY-method}

\title{Methods for Function \code{write.spectronaut} in Package \pkg{specL}}

\description{

Methods for function \code{write.spectronaut} in package \pkg{specL} ~~
writes specL objects to a file in a format which can be read 
by the 'Spectronaut' software.
additional arguments are

  \describe{

      \item{file}{A file name. default is \code{file='spec.txt'}. }

  }

}

\section{Methods}{

\describe{

\item{\code{signature(x = "specL")}}{
    Prints specL object data to to a file.
}

\item{\code{signature(x = "specLSet")}}{
    Prints specL object data to to a file.
}

}}

\keyword{methods}
