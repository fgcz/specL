\name{cdsw}

\alias{cdsw}
\alias{generateDynamicSwathWindow}


\title{Generate Dynamic SWATH Window}

\description{

    This function computes dynamic SWATH windows (\code{cdsw}) for a given
    vector of numeric values, an \code{psmSet} class, or an \code{specLSet} class.
    The input R data object can be generated using the 
    \code{read.bibliospec} function.
    
}

\usage{
    cdsw(x, n=20, overlap=1, ...)
}

\arguments{

  \item{x}{Numeric vector or psmSet class.}

  \item{n}{Number of desired SWATH windows. Default is set to 20.}

  \item{overlap}{Overlap of SWATH windows. The default is 1 Dalton.}

  \item{...}{pass arguments to \code{hist} function.}
}

\details{

  The function determines the SWATH windows using the quantile function.
  
}

\value{

  The output is \code{data.frame} having the attributes
  from, to, mid, width, and counts.
  In the ideal output all bins should contain the same number of numeric values.
  This requirement is violated because the window borders are rounded with no 
  digit after the comma.
  
}

\seealso{ 

    The S3 class definition:
  \code{showClass("psmSet")}

}

\author{

   Christian Panse, Christian Trachsel, 2015, 2016

}

\examples{

  # do not plot histogram
  cdsw(peptideStd, plot = FALSE,  overlap = 0)
  
  # plot hist
  cdsw(peptideStd, freq = TRUE)
  
  # pre-filtering
  ## cdsw(x=q1[350 <= q1 & q1 <= 1250], n=20, overlap = 0, freq=TRUE)
}
