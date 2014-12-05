#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/swath.R $
# $Id: swath.R 6657 2014-09-10 12:21:01Z cpanse $
# $Date: 2014-09-10 14:21:01 +0200 (Wed, 10 Sep 2014) $

# this function is for normalizing the rt on data
# for building the model data.fit  is used

.normalize_rt <- function(data, data.fit, iRT, plot=FALSE){

    message("normalizing RT ...")

    rt <- unlist(lapply(data, function(x){x$rt}))
    rt.fit <- unlist(lapply(data.fit, function(x){x$rt}))

    peptide <- as.character(unlist(lapply(data, function(x){x$peptideSequence})))
    peptide.fit <- as.character(unlist(lapply(data.fit, function(x){x$peptideSequence})))

    fileName <-  as.factor(unlist(lapply(data, function(x){x$fileName})))
    fileName.fit <-  as.factor(unlist(lapply(data.fit, function(x){x$fileName})))

    df <- data.frame(rt, peptide, fileName)
    df.fit <- data.frame(rt=rt.fit, peptide=peptide.fit, fileName=fileName.fit)

    data<-aggregate(df$rt, by=list(df$peptide, df$fileName), FUN=mean)
    data.fit<-aggregate(df.fit$rt, by=list(df.fit$peptide, df.fit$fileName), FUN=mean)

    names(data)<-c('peptide', 'fileName', 'aggregateInputRT')
    names(data.fit)<-c('peptide', 'fileName', 'aggregateInputRT')

    m <- merge(iRT, data.fit, by.x='peptide', by.y='peptide')

    if (nrow(m) < 3){
        message('not enough iRT peptides found for building the model.')
        message('=> no iRT regression applied, using orgiginal rt instead!')
        return(rt)
    }

    for (i in sort(unique(m$fileName))){message(paste("found", nrow(m[m$fileName == i,]), "iRT peptide(s) in", i)) }

    # build the model
    # we can do this only if we have found iRT peptides!

    nFileName <- length(unique (as.factor(fileName)))
    message("building model ...")

    if ( nFileName == 1){
	    message('model with only one file.')
    	fit <- lm(formula = rt ~ aggregateInputRT, data=m)

    } else if (nFileName > 1){
    	fit <- lm(formula = rt ~ aggregateInputRT * fileName, data=m)
    } else {
	    stop("problem in .normalize_rt.")
	}

    #fit$call

    # apply the model to my data 
    fileName[!fileName %in% m$fileName] <- NA
    rt.predicted <- predict(fit, data.frame(aggregateInputRT=rt, fileName=fileName))

    return(rt.predicted)
}

.defaultSwathFragmentIon <- function (b, y) {
        Hydrogen <- 1.007825
        Oxygen <- 15.994915
        Nitrogen <- 14.003074
            
        #bn_ <- (b + (n-1) * Hydrogen) / n 

        b1_ <- (b )
        y1_ <- (y ) 

        b2_ <- (b + Hydrogen) / 2
        y2_ <- (y + Hydrogen) / 2 

        b3_ <- (b + 2 * Hydrogen) / 3
        y3_ <- (y + 2 * Hydrogen) / 3

        return( cbind(b1_, y1_, b2_, y2_, b3_, y3_) )
}


genSwathIonLib <- function(data, 
    data.fit = data,
    mascotIonScoreCutOFF=20, 
    proteinIDPattern='', 
    max.mZ.Da.error = 0.1, 
    ignoreMascotIonScore = TRUE, 
    topN = 10,
    fragmentIonMzRange = c(200, 2000),
    fragmentIonRange = c(2,100), 
    fragmentIonFUN = .defaultSwathFragmentIon, 
    iRT = specL::iRTpeptides,
    AminoAcids = protViz::AA,
    file = NULL){ 

    if (fragmentIonRange[1] < 2){
        fragmentIonRange = c(2,100)
        warning("min fragmentIonRange should be at least set to 2. reset fragmentIonRange = c(2,100).")
    }

.genSwathIonLibSpecL <- function(x, fi, findNN.idx, mZ.error_, rt){
        m <- length(2 * nchar(x$peptideSequence))
        q1 <- x$pepmass
        q3 <- x$mZ[findNN.idx]
        fi.unlist<-sort(unlist(fi))
        q3.in_silico <- fi.unlist[ findNN_(q3, fi.unlist) ]

        if (sum(is.na(q3)) > 0){
            stop("ERROR")
        }

        irt <- round(rep(rt, m), 2)
        decoy <- rep(0, m)
        prec_z <- rep(x$charge, m)

        frg_type <- gsub("(.*)([0-9]+)_([0-9]+)", "\\1",  names(unlist(fi)))
        frg_z <- gsub("(.*)([0-9]+)_([0-9]+)", "\\2",  names(unlist(fi)))
        frg_nr  <- rep(1:nrow(fi), ncol(fi))  

        group_id <- rep(paste(x$peptideSequence, ".", x$charge,";", x$pepmass, sep=''), 1)

        # exspect a modification information, e.g., AAAMASATTM[+16.0]LTTK
        if ("peptideModSeq" %in% names(x)){
            if (nchar(x$peptideModSeq) > nchar(x$peptideSequenc)){
                group_id <- rep(paste(x$peptideModSeq, ".", x$charge,";", x$pepmass, sep=''), 1)
            }
        } else{
           # warning("x$peptideModSeq does not exists!")
        }

        intensity <- 100 * round(x$intensity[findNN.idx] / max(x$intensity[findNN.idx], na.rm=TRUE), 2)

        peptide_sequence <- rep(x$peptideSequence, 1)

        peptideModSeq <- x$peptideModSeq

        massErrorFilter <- ( (mZ.error_ < max.mZ.Da.error) & (fragmentIonMzRange[1] < q3 & q3 < fragmentIonMzRange[2]) )

        
        intensity.idx <- rev(order(intensity[massErrorFilter]))

        if (length(intensity.idx) < topN){
            topN <- length(intensity.idx)
        }
      
        idx <- intensity.idx[1:topN]

        res <- specL(group_id=group_id, 
                      peptide_sequence=peptide_sequence, 
                      proteinInformation=x$proteinInformation,
                      q1=x$pepmass, 
                      q3=as.numeric(q3[massErrorFilter])[idx], 
                      q3.in_silico=as.numeric(q3.in_silico[massErrorFilter])[idx], 
                      decoy=as.character(decoy)[idx], 
                      prec_z=x$charge, 
                      frg_type=frg_type[massErrorFilter][idx], 
                      frg_nr=as.numeric(frg_nr[massErrorFilter])[idx], 
                      frg_z=as.numeric(frg_z[massErrorFilter])[idx], 
                      relativeFragmentIntensity=as.numeric(intensity[massErrorFilter])[idx], 
                      irt=as.numeric(irt), 
                      peptideModSeq=as.character(peptideModSeq), 
                      mZ.error=as.numeric(mZ.error_[massErrorFilter])[idx], 
                      filename=x$fileName)
        
} # .genSwathIonLibSpecL

    x.peptideSeq<-unlist(lapply(data, function(x){x$peptideSequence}))

    x.varModMass<-(lapply(data, function(x){
        if (length(x$varModification) != nchar(x$peptideSequence)){
            rep(0, nchar(x$peptideSequence))
        }else{
            x$varModification
            }
        }))

    x.AAmass <- protViz::aa2mass(x.peptideSeq, AminoAcids$Monoisotopic, AminoAcids$letter1)

    x.AAmodifiedMass <- mapply(function(x,y){x+y}, x.varModMass, x.AAmass, SIMPLIFY = FALSE)

    x.rt <- unlist(lapply(data, function(x){x$rt}))

    if (length(iRT) > 1 & length(data.fit) > 1){
            x.rt <- .normalize_rt(data, data.fit, iRT, plot=FALSE)
    }

    message("generating ion library ...")
    # determine b and y fragment ions while considering the var mods
    fi<-lapply(x.AAmodifiedMass, function(x){fragmentIon(x, fragmentIonFUN)[[1]]})
    fragmentIonTyp = names(fi[[1]]) 
    
    # find NN peak
    findNN.idx<-mapply(function(x.fi, y.data){ findNN_(unlist(x.fi), y.data$mZ) }, 
        fi, data, 
        SIMPLIFY = FALSE)

    # determine mZ error
    mZ.error<-mapply(function(x, y.findNN.idx, z){
                abs(x$mZ[y.findNN.idx] - unlist(z))
            }, data, findNN.idx, fi, SIMPLIFY = FALSE)

    message("start generating specLSet object ..." )
    time.start <- Sys.time(); 

    # prepare table for output
    if (require(BiocParallel)){
        message("using BiocParallel::bpmapply( ..." )
        output <- BiocParallel::bpmapply (.genSwathIonLibSpecL, 
            data, fi, findNN.idx, mZ.error, x.rt, 
            SIMPLIFY = FALSE)
    }else{
        output <- mapply (.genSwathIonLibSpecL, 
            data, fi, findNN.idx, mZ.error, x.rt, 
            SIMPLIFY = FALSE)
    }

    time.end <- Sys.time();
    message(paste("time taken: ",  difftime(time.end, time.start, units='secs'), "secs"))

   
    return(specLSet(ionlibrary=output, 
        rt.normalized=unlist(x.rt), 
        rt.input=unlist(lapply(data, function(x){x$rt}))))
}

