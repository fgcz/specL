\name{NEWS}

\title{News for Package 'specL'}
\newcommand{\ghpr}{\href{https://github.com/protViz/specL/pull/#1}{##1}}
\newcommand{\ghit}{\href{https://github.com/protViz/specL/issues/#1}{##1}}

\section{CHANGES UNTIL VERSION 1.13.02 (2017-11-20)}{
\itemize{
	\item Eliminated C++-11 and Rcpp linking. Passing unit tests.
	\item use \code{lower_bound_} function exported by \url(https://cran.r-project.org/package=protViz)
}
}
\section{CHANGES UNTIL VERSION 1.9.15 (2017-04-21)}{
\itemize{
  \item added \CRANpkg{prozor} \code{(>= 0.2.2)} to the Suggest list.
  \item added more specific R package version numbers in DESCRIPTION file.
  \item in \code{plot.specLSet} (normalized RT versus RT) use \code{pch=16}
  and color with parameter \code{alpha=0.1}. 
  \item fixed issue \ghit{22} by including the iRTs in the ionlibrary; \code{LIB <- genSwathIonLib(data=peptideStd, data.fit=peptideStd.redundant); LIB@input.parameter$iRTpeptides}.
  \item fixed issue \ghit{19}.
  \item removed \code{par} command in specLset plot function.
  \item added \code{vignettes/report.Rmd} file, 
  see also \url{http://bioconductor.org/packages/devel/bioc/vignettes/specL/inst/doc/report.html}.
  }
}


\section{CHANGES IN VERSION 1.7.4 (2016-05-19)}{
  \itemize{
  \item USER VISIBLE CHANGES
    \itemize{
    \item added to \code{specLSet} summary "which std peptides (iRTs) where found in which raw files"
    \item one plot per raw file in plot methode of specLSet object
    }
  }
}

\section{CHANGES IN VERSION 1.7.1 (2016-05-13)}{
  \itemize{
  \item USER VISIBLE CHANGES
    \itemize{
    \item replaced NEWS by NEWS.Rd file
    \item modified specL object /replace decoy by score attribute \ghit{1},\ghit{4}
    \item specL on BioC 3.3
    }
  }
}


\section{CHANGES IN VERSION 1.5.10-13}{
  \itemize{
  \item USER UNVISIBLE CHANGES
    \itemize{
    \item added specLSet class support for cdsw methode
    \item changed Rmd5 vignette style
    \item added cdsw test case
    \item intro new vignette for cdsw method \url{http://bioconductor.org/packages/release/bioc/vignettes/specL/inst/doc/cdsw.html}
    }
  }
}

\section{CHANGES IN VERSION 1.5.9}{
    \itemize{
\item USER UNVISIBLE CHANGES
    \itemize{
  \item added test case for read.biliospec
    }}}
  
\section{CHANGES IN VERSION 1.5.5 - 1.5.8}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
  \item added RT prediction vignette file
  }}}

\section{CHANGES IN VERSION 1.5.4}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
  \item changed NAMESPACES and read.bibliospec docu to avoid warnings in 3.3 check
}}}

\section{CHANGES IN VERSION 1.5.3}{
\itemize{
\item USER UNVISIBLE CHANGES
  \itemize{
  
  \item added sqlite files for peptideStd RData \ghit{13}
  }}
}

\section{CHANGES IN VERSION 1.5.2}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
  \item find all signals having two or more in-silico fragment ions. \ghit{8}
  \item keep only the nearest fragment ion; if there are more take the first in line
}}}

\section{CHANGES IN VERSION 1.3.7}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
  \item getProteinPeptideTable - added
}}}

\section{CHANGES IN VERSION 1.3.5}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
  \item read.bibliospec - bugfixes
}}}

\section{CHANGES IN VERSION 1.3.4}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
  \item added Witold Wolski as maintainer
}
\item USER UNVISIBLE CHANGES
\itemize{
  \item read.bibliospec - replaced old code (for loop) by using mcmapply
  \item added time meassurements to read.bibliospec
}}}

\section{CHANGES IN VERSION 1.3.3}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
  \item plot::specLSet draws alpha circles iff plot(..., art=TRUE)
}

\item USER UNVISIBLE CHANGES
\itemize{
  \item .mascot2psmSet buxfix 
  \item renamed column  name in spectronaut outpu from irt to irt_or_rt
}}}

\section{CHANGES IN VERSION 1.3.2}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
  \item added ssrc (Sequence Specific Retention Calculator) function
  \item added a CITATION file
}}}

\section{CHANGES IN VERSION 1.3.1}{
\itemize{
\item USER VISIBLE CHANGES
  \itemize{
  \item added fucntion cdsw
}

\item USER UNVISIBLE CHANGES
\itemize{
  \item modified unit test for genSwathIonLib
}}}


\section{CHANGES IN VERSION 1.1.17}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
    \item removed file argument in genSwathIonLib function
}}}

\section{CHANGES IN VERSION 1.1.16}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
    \item added unit test for genSwathIonLib
}}}

\section{CHANGES IN VERSION 1.1.15}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
    \item added circle plots to specLSet plot method
    \item added breaks argument in genSwathIonLib methode
}}}


\section{CHANGES IN VERSION 1.1.14}{
\itemize{
\item USER UNVISIBLE CHANGES
\itemize{
  \item LinkedTo Rcpp; added C++ STL lower bound function 
  which is reqired for determining overlapping q1 and q3 SWATH windows
}}}

\section{CHANGES IN VERSION 1.1.13}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
  \item fixed man pages 
}}}

\section{CHANGES IN VERSION 1.1.12}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
  \item impoved package vignette
}}}

\section{CHANGES IN VERSION 1.1.11}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item modified default parameters of genSwathIonLib
    \item add content to vignette 
}}}

\section{CHANGES IN VERSION 1.1.10}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item added generate.consensus
}}}

\section{CHANGES IN VERSION 1.1.9}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item new features in specLSet 'summary' plot
}
\item USER UNVISIBLE CHANGES
\itemize{
    \item refactored merge.specLSet; merge by group_id 
    \item added unit test for merge.specLSet
}}}


\section{CHANGES IN VERSION 1.1.8}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item renamed annotateProteinID to annotate.protein_id
    \item added graphics on plot.specLSet method
}
\item USER UNVISIBLE CHANGES
\itemize{
    \item refactored merge
    }}}

\section{CHANGES IN VERSION 1.1.7}{
\itemize{

\item USER VISIBLE CHANGES
\itemize{
    \item introduce peakplot for bibliospec object
    \item introduce LCMS map for bibliospec object
    \item vignette cosmetics
}}}

\section{CHANGES IN VERSION 1.1.6}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item introduce specL_bibliospec summary method
}}}

\section{CHANGES IN VERSION 1.1.5}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item specLSet merge function 

    \item work on specLSet summary method
}}}

\section{CHANGES IN VERSION 1.1.4}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item summary method of specLSet class 
}
\item USER UNVISIBLE CHANGES
\itemize{
    \item unit test for data containing no iRT peptides
}}}

\section{CHANGES IN VERSION 1.1.3}{
\itemize{

\item USER VISIBLE CHANGES
\itemize{
    \item renamed write.Spectronaut to write.spectronaut

    \item write.spectronaut writes filename 

    \item added benchmark section in package vignette
}}}

\section{CHANGES IN VERSION 1.1.2}{

\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item uses modSeq in group_id iff existing
}}}

\section{CHANGES IN VERSION 1.1.1}{
\itemize{

\item USER VISIBLE CHANGES
\itemize{
    \item streamline modsequence, e.g., AAAMASATTM[+16.0]LTTK 
    for compatibility with peakView V2.0
}}}

\section{CHANGES IN VERSION 0.99.23}{
\itemize{

\item USER VISIBLE CHANGES
\itemize{
    \item added methods for specLSet class: 
        ionlibrary, rt.input, rt.normalized

    \item fixed Sys.time() units in message.
}
\item USER UNVISIBLE CHANGES
\itemize{
    \item genSwathIonLib using bpmapply
}}}


\section{CHANGES IN VERSION 0.99.22}{
\itemize{
\item USER VISIBLE CHANGES
\itemize{
    \item specLSet plot method
}}}

\section{CHANGES IN VERSION 0.99.21}{
\itemize{

\item USER VISIBLE CHANGES
\itemize{
    \item specLSet class

    \item replace print by show and write.Spectronaut method in specL and specLSet classes
    }}}
