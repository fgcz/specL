\docType{data}
\name{peptideStd}

\alias{peptideStd.redundant}
\alias{peptideStd}

\title{Peptide standard}

\format{contains a list of peptide spectrum assignments.}

\description{

    This dataset is a list of a peptide spectrum matches  
    (protein identification result) from two standard runs. 


}

\details{

    These standard runs (LCMS experiments) are routinely run on 
    well maintained instruments. In this case a standard run 
    consits of a digest of the FETUIN_BOVINE protein 
    (400 amol) and iRT peptides.

}

\author{

       Christian Panse, Christian Trachsel and Jonas Grossmann 2014

}

\examples{

    peakplot(peptideStd[[40]]$peptideSequence, peptideStd[[40]])

}


