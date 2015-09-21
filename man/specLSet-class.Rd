\name{specLSet-class}
\Rdversion{1.1}
\docType{class}

\alias{specLSet-class}

\alias{show,specLSet-method}
\alias{summary,specLSet-method}

\alias{plot}
\alias{plot,specLSet-method}

\alias{merge.specLSet}
\alias{merge.specLSet,specLSet-method}

\alias{write.spectronaut,specLSet-method}

\alias{ionlibrary}
\alias{ionlibrary,specLSet-method}

\alias{rt.input}
\alias{rt.input,specLSet-method}

\alias{rt.normalized}
\alias{rt.normalized,specLSet-method}

\alias{generate.consensus}
\alias{generate.consensus,specLSet-method}

\alias{ionlibrary}
\alias{ionlibrary,specLSet-method}

\alias{getProteinPeptideTable}
\alias{getProteinPeptideTable,specLSet-method}

\title{Class \code{"specLSet"}}

\description{

    This class is used to store, show, and write 
    generated results of the package.

}

\section{Objects from the Class}{

    Objects can be created by calls of the form \code{new("specLSet", ...)}.

}

\section{Slots}{

  \describe{

    \item{\code{input.parameter}:}{A list of parameter values. }
    
    \item{\code{ionlibrary}:}{A list of \code{"specL"} objects.}

    \item{\code{rt.normalized}:}{A numeric vector of  normalized retention time values.}

    \item{\code{rt.input}:}{A numeric vector of retention time values. }

  }

}

\section{Methods}{

  \describe{
    
    \item{show}{\code{signature(x = "specLSet")}: shows the object content. }

    \item{summary}{\code{signature(x = "specLSet")}: print summary of object content. }

    \item{plot}{\code{signature(x = "specLSet")}: plots normalized versus input rt. }
    
    \item{write.spectronaut}{\code{signature(x = "specLSet")}: writes object to ASCII file. }
    
    \item{generate.consensus}{\code{signature(x = "specLSet")}: generates consensus specLSet object by combining all specL objects having a redundant group_id which is defined by \code{paste(x@q1.in_silico, x@peptideModSeq, sep='_')}. }

    \item{merge.specLSet}{\code{signature(x = "specLSet")}: merges two specLSet objects.}

    \item{ionlibrary}{\code{signature(x = "specLSet")}: returns a list of \code{specL} objects. }

    \item{rt.input}{\code{signature(x = "specLSet")}: returns a numeric vector of input rt values. }

    \item{rt.normalized}{\code{signature(x = "specLSet")}: returns a numeric vector of normalized rt values. }
    
    \item{getProteinPeptideTable}{\code{signature(x = "specLSet")}: returns table of peptide protein mappings. }
	
    }

}

\author{

     Christian Panse 2014

}

\note{

    No notes yet.
}


\seealso{

     \code{\link{genSwathIonLib}}

}

\examples{

    showClass("specL")
    showClass("specLSet")

}
\keyword{classes}
