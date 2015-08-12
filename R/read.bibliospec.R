#R

# $HeadURL$
# $Id$
# $Date$

# TODO(cp): write a summery methide for read.bibliospec
# o list the number of files
# o list the number of specs

# usage:
# mZlist <- mcmapply(.decompose_peakMZ, data$peakMZ, data$numPeaks)
# maps base64 blob to human readable numbers
.decompose_peakMZ <- function(peakMZBlob, numPeaks) {
        mZ <- try(readBin(memDecompress(as.raw(peakMZBlob),'g'), double(), numPeaks), TRUE)
        if (!is.numeric(mZ)){
            mZ<-try(readBin(as.raw(peakMZBlob), double(),numPeaks), FALSE)
        }
        return(mZ)
}

# usage:
# intentityList <- mcmapply(.decompose_intensities, data$peakMZ, data$numPeaks)
# maps base64 blob to human readable numbers
.decompose_intensity <- function(peakIntensityBlob, numPeaks) {
        intensity <- try(readBin(memDecompress(as.raw(peakIntensityBlob),'g'), numeric(), n=numPeaks, size = 4), TRUE)
       
        if (!is.numeric(intensity) || length(intensity) != numPeaks){
            intensity <- try(readBin(as.raw(peakIntensityBlob), 
                numeric(), n=numPeaks, size = 4), FALSE)
        }  
        return (intensity)
}

.convert_blib2psm_parallel <- function(data, ncores=1){
  mZ_list <- mcmapply(specL:::.decompose_peakMZ, data$peakMZ, data$numPeaks, mc.cores=ncores, mc.preschedule=TRUE)
  message(paste("decomposed mZ values"))
  
  intentity_list <- mcmapply(specL:::.decompose_intensity, data$peakIntensity, data$numPeaks, mc.cores=ncores, mc.preschedule=TRUE)
  message(paste("decomposed intensity values"))
  
  res <- mcmapply(function(peaks, mZ, intensity, peptideSequence, peptideModSeq, 
                         charge, pepmass, fileName, rt, score, mc.cores=ncores, mc.preschedule=TRUE){
    psm <-list(peaks=peaks, 
               mZ=mZ, 
               intensity=intensity,
               peptideSequence=peptideSequence,
               peptideModSeq=peptideModSeq,
               charge=charge, 
               pepmass=pepmass,
               fileName = fileName,
               proteinInformation='',
               rt=rt,
               varModification=rep(0.0, nchar(peptideSequence)),
               mascotScore = -10 * log((1E-6 + score)) / log(10)
               );
    
    class(psm) <- "psm"
    return(psm)
    },
    data$numPeaks,
    mZ_list, 
    intentity_list,
    data$peptideSeq,
    data$peptideModSeq,
    data$precursorCharge,
    data$precursorMZ,
    data$fileName,
    data$retentionTime,
    data$score,
    SIMPLIFY = FALSE, mc.cores=ncores, mc.preschedule=TRUE)
  
  class(res) <- 'psmSet'
  return(res)
}

.convert_blib2psm <- function(data){
    res <- list()

      
     for (i in 1:nrow(data)){
        x <- data[i, ]

        mZ <- try(readBin(memDecompress(as.raw(x$peakMZ[[1]]),'g'), double(), x$numPeaks), TRUE)

        if (!is.numeric(mZ)){
            mZ<-try(readBin(as.raw(x$peakMZ[[1]]), double(),x$numPeaks), FALSE)
        }

        intensity <- try(readBin(memDecompress(as.raw(x$peakIntensity[[1]]),'g'), 
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
        class(res[[i]]) = "psm"
    }
	return(res)
}



read.bibliospec <- function(file,ncores=NULL){
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

    data <- DBI::fetch(SQLQuery0, n = -1)

    if (msg<-dbGetException(con)$errorNum != 0){
        stop(msg$errorMsg)
    }

    SQLQuery1 <- dbSendQuery(con, 
        statement = "SELECT RefSpectraID, position, mass FROM Modifications;")

    if (msg<-dbGetException(con)$errorNum != 0){
        stop(msg$errorMsg)
    }

    data.modifications <- DBI::fetch(SQLQuery1, n = -1)

    if (msg<-dbGetException(con)$errorNum != 0){
        stop(msg$errorMsg)
    }


    message(paste("fetched", nrow(data), "rows."))
    
    res<-list()
    if (require(parallel)){
      if(is.null(ncores)){
        ncores <- max(1,parallel::detectCores()/2)
      }
      message(paste("start converting blib blobs to psm using", ncores, "cores ..."))
      time.start <- Sys.time(); 
      res <- .convert_blib2psm_parallel(data, ncores)
      time.end <- Sys.time();
      message(paste("time taken: ", round(difftime(time.end, time.start, units='mins'),2),  "minutes"))
      
    }else{
      # TODO(cp): replace it after testing if the result of the methode used in the
      # block above is the same
      res <- .convert_blib2psm(data)
    }

    message(paste("assigning", nrow(data.modifications), "modifications ..."))
    for (i in 1:nrow(data.modifications)){
        res[[data.modifications$RefSpectraID[i]]]$varModification[data.modifications$position[i]] <- data.modifications$mass[i]
    }

    class(res)='psmSet'

    dbDisconnect(con)
    return(res)
}
#s<-read.bibliospec("/scratch/specL_revisions_201412/p1000_testBelowFour.redundant.blib")

summary.psmSet <- function (object, ...){
  
    cat("Summary of a \"psmSet\" object.")

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

plot.psmSet <- function (x, iRTpeptides=specL::iRTpeptides, ...){
  
  rt <- unlist(lapply(x, function(x){x$rt}))
  pepmass <- unlist(lapply(x, function(xx){xx$pepmass}))
  
  peptide <- unlist(lapply(x, function(xx){xx$peptideSequence}))
  idx.iRT <- which(peptide %in% iRTpeptides$peptide)
  
  charge <- unlist(lapply(x, function(xx){xx$charge}))
  filename <- as.numeric(as.factor(unlist(lapply(x, function(xx){xx$fileName}))))
  
  plot(pepmass ~ rt, 
       pch=filename, 
       col=charge, 
       main='LCMS map',
       xlab='retention time',
       ylab='peptide mass', ...)
  
  cc<-sort(unique(charge))
  legend('topleft', "iRT peptides", pch=22)
  legend('left', paste(cc,'+',sep=''), col=cc, pch=22, title='charge')
  
  text(rt,pepmass, 1:length(rt),pos=3,col=charge,cex=0.5)
  
  fn<-unique(unlist(lapply(x, function(xx){xx$fileName})))
  n<-nchar(fn)
  legend('bottomright', 
         substr(fn, n - 25, n), 
         pch=unique(filename),
         cex=1.0, title='input file names')
  
  points(rt[idx.iRT], pepmass[idx.iRT], pch=22, cex=2)
  
}

plot.psm <- function (x, ...){
  
  AAmass <- protViz::aa2mass(x$peptideSequence)[[1]]#, protViz::AA$Monoisotopic, protViz::AA$letter1)
  
  AAmodifiedMass <- AAmass + x$varModification

  fi <- protViz::fragmentIon(AAmodifiedMass, FUN=.defaultSwathFragmentIon)[[1]]
  
  spec <- list(mZ=x$mZ, intensity=x$intensity)
  
  return(protViz::peakplot(peptideSequence=x$peptideSequence, spec=spec, fi=fi, ...))
}

.mascot2psmSet <-  function(dat, mod, mascotScoreCutOff=40){
    res <- lapply(dat, function(x){
      x$MonoisotopicAAmass <- protViz::aa2mass(x$peptideSequence)[[1]]#, protViz::AA$Monoisotopic, protViz::AA$letter1)
      
      modString <- as.numeric(strsplit(x$modification, '')[[1]])
      modIdx <- which(modString > 0.0) - 1
      modString.length <- length(modString)
      
      x$varModification <- mod[modString [c(-1, -modString.length)] + 1 ] 
      if (length(modIdx) > 0){
        warning("modified varModification caused.")
        x$varModification[modIdx] <- x$varModification[modIdx] - x$MonoisotopicAAmass[modIdx]
      }
        rt<-x$rtinseconds
      x<-c(x, rt=rt, fileName="mascot")
      
      class(x) <- "psm"
      return(x)
    })
    
    # filter
    res<-res[which(unlist(lapply(dat, function(x){
      x$mascotScore > mascotScoreCutOff && length(x$mZ)>10}))) ]
    class(res) <- "psmSet"
    return(res)
  }
