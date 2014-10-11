#R

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
            charge=x$precursorCharge, 
            pepmass=x$precursorMZ,
            fileName = x$fileName,
            proteinInformation='',
            rt=x$retentionTime,
            varModification=rep(0.0, nchar(x$peptideSeq)),
            mascotScore = -10 * log((1E-6 + x$score)) / log(10))
    }

    message(paste("assigning", nrow(data.modifications), "modifications ..."))
    for (i in 1:nrow(data.modifications)){
        res[[data.modifications$RefSpectraID[i]]]$varModification[data.modifications$position[i]] <- data.modifications$mass[i]
    }

    dbDisconnect(con)
    return(res)
}
