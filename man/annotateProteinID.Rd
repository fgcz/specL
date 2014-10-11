\name{annotateProteinID}

\alias{annotateProteinID}

\title{annotateProteinID}

\description{

This function assigns the protein identifier
for a list of tandem mass specs having
a peptide sequence assigned.

}

\usage{

annotateProteinID(data, file = NULL, fasta = read.fasta(file = file, 
         as.string = TRUE, seqtype = "AA"), digestPattern = "(([RK])|(^)|(^M))") 

}

\arguments{

  \item{data}{list of records containing mZ and peptide sequences.}

  \item{file}{file name of a FASTA file.}

  \item{fasta}{a fasta object as reaturned by the \code{seqinr::read.fasta(...)} 
  methode.}

  \item{digestPattern}{a regex pattern which can be used by the \code{grep} 
  command. the default regex pattern assumes a tryptic digest.}

}

\details{

The protein sequences a read by the \code{read.fasta} function
of the \code{seqinr} package. The protein identifier is written 
to the protein \code{proteinInformation} variable.

If the function is called on a multi core architecture it uses \code{mclapply}.

It is recommended to load the FASTA file prior to running 
\code{annotateProteinID} using \code{myFASTA <- read.fasta(file = file, 
as.string = TRUE, seqtype = "AA")} instead of providing the FASTA file name to 
the function.
}

\seealso{

  \code{?read.fasta} of the \code{seqinr} package.
  
  \url{http://www.uniprot.org/help/fasta-headers}
  
}


\author{

   Jonas Grossmann and Christian Panse, 2014
}

\value{

    it returns a list object.

}

\examples{

    # annotateProteinID
    
    # our Fasta sequence
      irtFASTAseq <- paste(">zz|ZZ_FGCZCont0260|", 
      "iRT_Protein_with_AAAAK_spacers concatenated Biognosys\n",
      "LGGNEQVTRAAAAKGAGSSEPVTGLDAKAAAAKVEATFGVDESNAKAAAAKYILAGVENS",
      "KAAAAKTPVISGGPYEYRAAAAKTPVITGAPYEYRAAAAKDGLDAASYYAPVRAAAAKAD",
      "VTPADFSEWSKAAAAKGTFIIDPGGVIRAAAAKGTFIIDPAAVIRAAAAKLFLQFGAQGS",
      "PFLK\n")
      
    # be realistic, do it from file
      Tfile <- file();  cat(irtFASTAseq, file = Tfile);
      
    #use read.fasta from seqinr
      fasta.irtFASTAseq <-read.fasta(Tfile, as.string=TRUE, seqtype="AA")
      close(Tfile)
    
    #annotate with proteinID 
    # -> here we find all psms from the one proteinID above
      peptideStd <- specL::annotateProteinID(peptideStd, 
      fasta=fasta.irtFASTAseq)
  
    #show indices for all PSMs where we have a proteinInformation
     which(unlist(lapply(peptideStd, 
      function(x){nchar(x$proteinInformation)>0})))

}
