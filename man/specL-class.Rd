\name{specL-class}
\Rdversion{1.1}
\docType{class}

\alias{specL-class}

\alias{plot,specL-method}

\alias{show,specL-method}

\alias{write.spectronaut,specL-method}

\title{Class \code{"specL"}}

\description{
    This class is used to store, print, and plot  the
       generated results of the package.

}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("specL", ...)}.
}

\section{Slots}{

  \describe{
    \item{\code{group_id}:}{Object of class \code{"character"} just an id }
    \item{\code{peptide_sequence}:}{Object of class \code{"character"} AA 
    sequence }
    
    \item{\code{proteinInformation}:}{Object of class \code{"character"} a 
    string contains the protein identifier. }
    
    \item{\code{q1}:}{Object of class \code{"numeric"} peptide weight m/Z as 
    measured by the MS device }
    
    \item{\code{q3}:}{Object of class \code{"numeric"} measured fragment ions. }
    \item{\code{q3.in_silico}:}{Object of class \code{"numeric"} in-silico 
    derived fragment ions. }
    
    \item{\code{decoy}:}{Object of class \code{"character"} is this  a decoy 
    hit? 1 or 0. }
    
    \item{\code{prec_z}:}{Object of class \code{"numeric"} pre-cursor charge. }
    \item{\code{frg_type}:}{Object of class \code{"character"} fragmenbt ion 
    type, e.g., b or y ion. }
    
    \item{\code{frg_nr}:}{Object of class \code{"numeric"} fragment ion number}
    \item{\code{frg_z}:}{Object of class \code{"numeric"} fragment ion charge. }
    \item{\code{relativeFragmentIntensity}:}{Object of class \code{"numeric"} 
    percentage base peaks of frament ions. }
    
    \item{\code{irt}:}{Object of class \code{"numeric"} independent retention 
    time in seconds. }
    
    \item{\code{peptideModSeq}:}{Object of class \code{"numeric"} a vector 
    contains the mass diff between AA and mod AA. }
    
    \item{\code{mZ.error}:}{Object of class \code{"numeric"} a string contains 
    the protein identifier. }

    \item{\code{filename}:}{Object of class \code{"character"} a string contains 
    the filename of the ions.}
  }
}

\section{Methods}{
  \describe{

    \item{plot}{\code{signature(x = "specL")}: plots the fragment ions of \code{specL} object.}

    \item{show}{\code{signature(x = "specL")}: shows the content of \code{specL} object.}

    \item{write.spectronaut}{\code{signature(x = "specL")}: writes the \code{specL} object to a ASCII file. }

	 }
}

%\references{             
%}

\author{

     Christian Panse 2014

}
\note{

    No notes yet.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

     \code{\link{genSwathIonLib}}

}
\examples{

    showClass("specL")

}
\keyword{classes}
