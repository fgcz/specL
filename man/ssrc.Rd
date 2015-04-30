\name{ssrc}

\alias{ssrc}
\alias{hydrophobicity}


\title{Sequence Specific Retention Calculator}

\description{

    This function returns as output one hydrophobicity value 
    for a given sequence of amino acids (tryptic peptide)
    which can be used to predict the retention times.
    The calculation is based on the method described in PMID:15238601.
    
}

\usage{

    ssrc(x, H=list())
}

\arguments{

  \item{x}{sequence of amino acids, e.g., \code{x="ELIVSK"}}

  \item{H}{A list of retention coefficients. The default is set to the values of PMID:15238601 table II column 2(Rc values).}

}

\details{

  t.b.d.
  
}

\seealso{ 

    Krokhin, O. V. et al. An improved model for prediction of retention times of tryptic peptides in ion
    pair reversed-phase HPLC: its application to protein peptide mapping by off-line HPLC-MALDI MS.
    Mol. Cell Proteomics 3, 908-919 (2004).
    \url{http://www.ncbi.nlm.nih.gov/pubmed/15238601}

    \url{http://hs2.proteome.ca/SSRCalc/SSRCalcX.html}

}

\author{

   Christian Panse, Christian Trachsel 2015

}

\examples{

  # example of table iv [PMID:15238601]
  lapply(c("SCHTAVGR", "SCHTGLGR", "EDLIAYLK"), ssrc)
  
  
    Rc <- list()
    Rc[['W']] <- 11.0; 
    Rc[['F']] <- 10.5; 
    Rc[['L']] <- 9.6; 
    Rc[['I']] <- 8.4; 
    Rc[['M']] <- 5.8; 
    Rc[['V']] <- 5.0; 
    Rc[['Y']] <- 4.0; 
    Rc[['A']] <- 0.8; 
    Rc[['T']] <- 0.4; 
    Rc[['P']] <- 0.2; 
    Rc[['E']] <- 0.0; 
    Rc[['D']] <- -0.5; 
    Rc[['C']] <- -0.8; 
    Rc[['S']] <- -0.8; 
    Rc[['Q']] <- -0.9; 
    Rc[['G']] <- -0.9; 
    Rc[['N']] <- -1.2; 
    Rc[['R']] <- -1.3; 
    Rc[['H']] <- -1.3; 
    Rc[['K']] <- -1.9; 

    RcNt <- lapply(Rc, function(x){sum(unlist(Rc))/20 - x})
    
    ssrc <- unlist(lapply(peptideStd, function(x){ssrc(x$peptideSequence)}))
    rt <- unlist(lapply(peptideStd, function(x){x$rt}))

    plot(ssrc, rt) 
    abline(ssrc.lm <- lm(rt ~ ssrc), col='red'); 
    legend("topleft", paste("spearman", round(cor(ssrc, rt, method='spearman'),2)))

}
