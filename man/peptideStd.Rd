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

    \dontrun{

        download.file("http://fgcz-data.uzh.ch/~cpanse/specL/data/peptideStd.blib", 
            destfile="peptideStd.blib")

        download.file("http://fgcz-data.uzh.ch/~cpanse/specL/data/peptideStd.redundant.blib", 
            destfile="peptideStd.redundant.blib")

        # checksum

        if (require(tools)){
            md5sum(c("peptideStd.blib", "peptideStd.redundant.blib")) == 
                c("3f231931e54efd6516d7aa302073b17f", 
                    "8bab829d9e99344136613a17c0374b90")
        }

        peptideStd <- read.bibliospec("peptideStd.blib")
        peptideStd.redundant <- read.bibliospec("peptideStd.redundant.blib")

    }

}


\seealso{

    \url{http://fgcz-data.uzh.ch/~cpanse/specL/data/}

}
