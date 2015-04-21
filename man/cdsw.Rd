\name{cdsw}

\alias{cdsw}
\alias{generateDynamicSwathWindow}


\title{Generate Dynamic SWATH Window}

\description{

    This function generates an ion library for SWATH analysis.
    It takes a R data object which contains a peak list.
    The R data object can be generated using the 
    \code{read.bibliospec} function.
    
}

\usage{
    cdsw(x, n=20, overlap=1)
}

\arguments{

  \item{x}{Numeric vector or psmSet class.}

  \item{n}{Number of desired SWATH windows.}

  \item{overlap}{Overlap of SWATH windows. The default is 1 Dalton}

}

\details{

  The function determines the SWATH windows using the quantile function.

}

\value{

  The output is table containing having the attributes
  from, to, mid, width, and counts.
  In the idea output all bins should contain the same number of numeric values.

}

\seealso{ 

    The S3 class definition:
  \code{showClass("psmSet")}

}

\author{

   Christian Panse, Christian Trachsel 2015

}

\examples{

  cdsw(peptideStd)
  
}
