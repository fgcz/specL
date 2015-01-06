#R

# $HeadURL$
# $Id$
# $Date$


# TODO(cp): write a summery methide for read.bibliospec
# o list the number of files
# o list the number of specs

read.bibliospec <- function(file){
    m <- dbDriver("SQLite", max.con=25)       
    con <- dbConnect(m , dbname=file, flags = SQLITE_RO)

    SQLQuery0 <- dbSendQuery(con, statement = paste(
        "SELECT numPeaks, peakMZ, peakIntensity, peptideSeq,",
        "precursorCharge, precursorMZ, retentionTime,",
        "peptideModSeq, score, SpectrumSourceFiles.fileName",
        "FROM SpectrumSourceFiles, RefSpectraPeaks, RefSpectra",
        "WHERE RefSpectra.id=RefSpectraPeaks.RefSpectraID",
        "and SpectrumSourceFiles.id = RefSpectra.fileID;", 
        sep=" "))


    if (msg<-dbGetException(con)$errorNum != 0){
        stop(msg$errorMsg)
    }

    data <- fetch(SQLQuery0, n = -1)

    if (msg<-dbGetException(con)$errorNum != 0){
        stop(msg$errorMsg)
    }

    SQLQuery1 <- dbSendQuery(con, 
        statement = "SELECT RefSpectraID, position, mass FROM Modifications;")

    if (msg<-dbGetException(con)$errorNum != 0){
        stop(msg$errorMsg)
    }

    data.modifications <- fetch(SQLQuery1, n = -1)

    if (msg<-dbGetException(con)$errorNum != 0){
        stop(msg$errorMsg)
    }


    message(paste("fetched", nrow(data), "rows."))
    res <- list()

     for (i in 1:nrow(data)){

        x <- data[i, ]

        mZ <- try(readBin(memDecompress(as.raw(x$peakMZ[[1]]),'g'), double(), x$numPeaks), TRUE)

        if (!is.numeric(mZ)){
            mZ<-try(readBin(as.raw(x$peakMZ[[1]]), double(),x$numPeaks), FALSE)
        }

        intensity<-try(readBin(memDecompress(as.raw(x$peakIntensity[[1]]),'g'), 
            numeric(), n=x$numPeaks, size = 4), TRUE)

        if (!is.numeric(intensity) || length(intensity) != length(mZ)){
            intensity <- try(readBin(as.raw(x$peakIntensity[[1]]), 
                numeric(), n=x$numPeaks, size = 4), FALSE)
        }

        if ( length(intensity) != length(mZ)  ){
            warning(" length(intensity) != length(mZ) ")
        }

        res[[i]] <- list(
            peaks=x$numPeaks,
            mZ=mZ, 
            intensity=intensity, 
            peptideSequence=x$peptideSeq,
            peptideModSeq=x$peptideModSeq,
            charge=x$precursorCharge, 
            pepmass=x$precursorMZ,
            fileName = x$fileName,
            proteinInformation='',
            rt=x$retentionTime,
            varModification=rep(0.0, nchar(x$peptideSeq)),
            mascotScore = -10 * log((1E-6 + x$score)) / log(10))
        class(res[[i]]) = "psm_bibliospec"
    }

    message(paste("assigning", nrow(data.modifications), "modifications ..."))
    for (i in 1:nrow(data.modifications)){
        res[[data.modifications$RefSpectraID[i]]]$varModification[data.modifications$position[i]] <- data.modifications$mass[i]
    }

    class(res)='psmSet_bibliospec'

    dbDisconnect(con)
    return(res)
}
#s<-read.bibliospec("/scratch/specL_revisions_201412/p1000_testBelowFour.redundant.blib")

summary.psmSet_bibliospec <- function (object, ...){
  
    cat("Summary of a \"specL_bibliospec\" object.")

    cat("\nNumber of precursor:\n\t")
    cat(length(object))

    cat("\nNumber of precursors in Filename(s)\n")
    t <- (table(unlist(lapply(object, function(x){x$fileName}))))
      
    n <- names(t)
    w <- getOption("width") / 2
    for (i in 1:length(t)){
      cat('\t')
      cat(substr(n[i], nchar(n[i])-w, nchar(n[i])))
      cat('\t')
      cat(t[i])
      cat('\n')
    }
    
    cat("Number of annotated precursor:\n\t")
    cat(sum(unlist(lapply(object, function(x){x$proteinInformation != ''}))))
    cat ("\n")
}

plot.psmSet_bibliospec <- function (object, ...){
  
  rt <- unlist(lapply(object, function(x){x$rt}))
  pepmass <- unlist(lapply(object, function(x){x$pepmass}))
  charge <- unlist(lapply(object, function(x){x$charge}))
  filename <- as.numeric(as.factor(unlist(lapply(object, function(x){x$fileName}))))
  
  plot(pepmass ~ rt, 
       pch=filename, 
       col=charge, 
       main='LCMS map',
       xlab='retention time',
       ylab='peptide mass', ...)
  
  cc<-sort(unique(charge))
  legend('topleft', paste(cc,'+',sep=''), col=cc, pch=22)
  
  text(rt,pepmass, 1:length(rt),pos=3,col=charge,cex=0.5)
  
  fn<-unique(unlist(lapply(object, function(x){x$fileName})))
  n<-nchar(fn)
  legend('bottomright', 
         substr(fn, n - 25, n), 
         pch=unique(filename),cex=1.0)
  
}

plot.psm_bibliospec <- function (object, ...){
  
  AAmass <- protViz::aa2mass(object$peptideSequence)[[1]]#, protViz::AA$Monoisotopic, protViz::AA$letter1)
  
  AAmodifiedMass <- AAmass + object$varModification

  fi <- protViz::fragmentIon(AAmodifiedMass, FUN=.defaultSwathFragmentIon)[[1]]
  
  spec <- list(mZ=object$mZ, intensity=object$intensity)
  
  return(protViz::peakplot(peptideSequence=object$peptideSequence, spec=spec, fi=fi, ...))
}
