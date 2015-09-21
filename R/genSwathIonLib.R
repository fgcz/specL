#R

# $HeadURL$
# $Id$
# $Date$



.plot_ptm_ion <- function(psm){
  # find all spec idxs having at least one modification
  psm.mod_idx <- which(lapply(psm, function(x){sum(x$varModification != 0.0)}) > 0)
 
  message(psm.mod_idx)
  # draw all spec with index psm.mod_idx
  lapply(psm.mod_idx, function(idx){
    pp <- plot(psm[[idx]], sub=psm[[idx]]$peptideModSeq)
    
    b_idxs <- which(psm[[idx]]$varModification != 0.0)

    nAA <- length(psm[[idx]]$varModification)
 
    y_idxs <- nAA - b_idxs + 1

    abline(v=pp$fragmentIon$b1_[b_idxs], col='purple')
    abline(v=pp$fragmentIon$y1_[y_idxs], col='cyan')
    
    abline(v=pp$fragmentIon$b2_[b_idxs], col='purple', lwd=2)
    abline(v=pp$fragmentIon$y2_[y_idxs], col='cyan', lwd=2)
  })
}


# .swath_plot_(peptideStd[[136]])
.swath_plot_ <- function(x){
  
  
  x.AAmass <- protViz::aa2mass(x$peptideSequence, protViz::AA$Monoisotopic, protViz::AA$letter1)
  
  x.AAmodifiedMass <- mapply( function(x, y){ x + y }, x$varModification, x.AAmass, SIMPLIFY = FALSE)
  # fi <- lapply( x.AAmodifiedMass, function(x){ fragmentIon(x, fragmentIonFUN)[[1]] } )
  
  fi <- protViz::fragmentIon(x.AAmodifiedMass)[[1]]
  idx <- which(x$varModification != 0.0)
  
  protViz::peakplot(peptideSequence=x$peptideSequence, spec=x)
  abline(v=fi$b[idx])
  abline(v=fi$y[idx], col='green')
  
  print(fi)
  print(idx)
}

.swath_q1q3_conflicts <- function(ionlib, breaks=seq(400,2000,by=25), overlap=1.0){
  
  
  lapply(ionlib, function(x){
    q3 <- x@q3
    q1 <- x@q1
  
    q1_idx <- .Call("lower_bound_", q1, breaks, PACKAGE = "specL")
    q3_idx <- .Call("lower_bound_", q3, breaks, PACKAGE = "specL")
  

    res <- cbind(q1, q1_idx, q3_idx,  
               lower=breaks[q3_idx], q3, upper=breaks[q3_idx+1], 
               check=(breaks[q3_idx] < rep(q3, length(q3_idx)) & rep(q3, length(q3_idx)) < breaks[q3_idx+1]), 
               conflict=(rep(q1_idx, length(q3_idx)) == q3_idx))
  
    #print(cbind(q1, q1_idx, breaks[q1_idx], breaks[q1_idx+1])) 
    # print(res)
  })

}


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

    if (nrow(m) < 1){
      message('no iRT peptides found for building the model.')
      message('=> no iRT regression applied, using orgiginal rt instead!')
      return(rt)
    }
    
    if (nrow(m) < 3){
        message('not enough iRT peptides found for building the model.')
        message('=> no iRT regression applied, using original rt instead!')
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
    fragmentIonMzRange = c(300, 1800),
    fragmentIonRange = c(5, 100), 
    fragmentIonFUN = .defaultSwathFragmentIon, 
    iRT = specL::iRTpeptides,
    AminoAcids = protViz::AA,
    breaks=NULL){ 

    # one transition is useless anyway
    if (fragmentIonRange[1] < 2){
        fragmentIonRange = c(2,100)
        warning("min fragmentIonRange should be at least set to 2. reset fragmentIonRange = c(2,100).")
    }

.genSwathIonLibSpecL <- function(x, fi, findNN.idx, mZ.error_, rt){
        m <- length(2 * nchar(x$peptideSequence))
        q1 <- x$pepmass
        q3 <- x$mZ[findNN.idx]
        fi.unlist <- sort(unlist(fi))
        q3.in_silico <- fi.unlist[ findNN_(q3, fi.unlist) ]
        q1.in_silico <- protViz::parentIonMass(x$peptideSequence) + sum(x$varModification)

        if (sum(is.na(q3)) > 0){
            stop("ERROR")
        }

        #message(paste("length of q3 ", length(q3)))

        irt <- round(rep(rt, m), 2)
        decoy <- rep(0, m)
        prec_z <- rep(x$charge, m)

        frg_type <- gsub("(.*)([0-9]+)_([0-9]+)", "\\1",  names(unlist(fi)))
        frg_z <- gsub("(.*)([0-9]+)_([0-9]+)", "\\2",  names(unlist(fi)))
        frg_nr  <- rep(1:nrow(fi), ncol(fi))  

        # group_id <- rep(paste(x$peptideSequence, ".", x$charge,";", x$pepmass, sep=''), 1)
        #group_id <- rep(paste(x$peptideSequence, ".", x$charge,";", round(q1.in_silico,3), sep=''), 1)
        group_id <- rep(paste(x$peptideSequence, ".", x$charge, sep=''), 1)
        
        # expect a modification information, e.g., AAAMASATTM[+16.0]LTTK
        if ("peptideModSeq" %in% names(x)){
            if (nchar(x$peptideModSeq) > nchar(x$peptideSequenc)){
                group_id <- rep(paste(x$peptideModSeq, ".", x$charge, sep=''), 1)
            }
        } else{
           # warning("x$peptideModSeq does not exists!")
        }

        intensity <- 100 * round(x$intensity[findNN.idx] / max(x$intensity[findNN.idx], na.rm=TRUE), 2)

        peptide_sequence <- rep(x$peptideSequence, 1)

        peptideModSeq <- x$peptideModSeq
        
        filter_mass_error <- (mZ.error_ < max.mZ.Da.error) & (fragmentIonMzRange[1] < q3 & q3 < fragmentIonMzRange[2])
        
        # TODO(cp): add the SWATH window filter here
        if (length(breaks) > 1){
          q1_idx <- .Call("lower_bound_", q1, breaks, PACKAGE = "specL")
          q3_idx <- .Call("lower_bound_", q3, breaks, PACKAGE = "specL")
          
          filter_swath_window <- (rep(q1_idx, length(q3_idx)) != q3_idx)
          if (length(filter_mass_error) != length(filter_swath_window)){
            warning("filter have different length!")
          }
          filter_mass_error <- filter_mass_error & filter_swath_window
        }
        
        #print (sum(filter_mass_error))
        intensity.idx <- rev(order(intensity[filter_mass_error]))
        
        if (length(intensity.idx) < topN){
            topN <- length(intensity.idx)
        }
      
        idx <- intensity.idx[1:topN]

        res <- specL(group_id=group_id, 
                      peptide_sequence=peptide_sequence, 
                      proteinInformation=x$proteinInformation,
                      q1=x$pepmass, 
                      q1.in_silico = q1.in_silico,
                      q3=as.numeric(q3[filter_mass_error])[idx], 
                      q3.in_silico=as.numeric(q3.in_silico[filter_mass_error])[idx], 
                      decoy=as.character(decoy)[idx], 
                      prec_z=x$charge, 
                      frg_type=frg_type[filter_mass_error][idx], 
                      frg_nr=as.numeric(frg_nr[filter_mass_error])[idx], 
                      frg_z=as.numeric(frg_z[filter_mass_error])[idx], 
                      relativeFragmentIntensity=as.numeric(intensity[filter_mass_error])[idx], 
                      irt=as.numeric(irt), 
                      peptideModSeq=as.character(peptideModSeq), 
                      mZ.error=as.numeric(mZ.error_[filter_mass_error])[idx], 
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

    x.AAmodifiedMass <- mapply(function(x, y){x + y}, x.varModMass, x.AAmass, SIMPLIFY = FALSE)

    x.rt <- unlist(lapply(data, function(x){x$rt}))

    if (length(iRT) > 1 & length(data.fit) > 1){
            x.rt <- .normalize_rt(data, data.fit, iRT, plot=FALSE)
    }

    message("generating ion library ...")
    # determine fragment ions while considering the var mods
    fi <- lapply(x.AAmodifiedMass, function(x){fragmentIon(x, fragmentIonFUN)[[1]]})

    fragmentIonTyp = names(fi[[1]]) 
    
    # find NN peak
    findNN.idx <- mapply(function(x.fi, y.data){ findNN_(unlist(x.fi), y.data$mZ) }, 
        fi, data, 
        SIMPLIFY = FALSE)

    # determine mZ error
    mZ.error <- mapply(function(x, y.findNN.idx, z){
                abs(x$mZ[y.findNN.idx] - unlist(z))
            }, data, findNN.idx, fi, SIMPLIFY = FALSE)

    message("start generating specLSet object ..." )
    time.start <- Sys.time(); 

    message(paste("length of findNN idx ", length(findNN.idx)))

    # prepare table for output
    if ( parallel::detectCores() > 1 & length(data) > 200){
        message("using BiocParallel::bpmapply( ..." )
        output <- parallel::mcmapply (.genSwathIonLibSpecL, 
            data, fi, findNN.idx, mZ.error, x.rt, 
            SIMPLIFY = FALSE)
    }else{
        output <- mapply (.genSwathIonLibSpecL, 
            data, fi, findNN.idx, mZ.error, x.rt, 
            SIMPLIFY = FALSE)
    }

    time.end <- Sys.time();
    message(paste("time taken: ",  difftime(time.end, time.start, units='secs'), "secs"))

    output <- output[which(unlist(lapply (output, function (x) {fragmentIonRange[1] <= length(x@q3) && length(x@q3) <= fragmentIonRange[2]})))]
   
    return(specLSet(ionlibrary=output, 
        input.parameter=list(mascotIonScoreCutOFF=mascotIonScoreCutOFF, 
                   proteinIDPattern=proteinIDPattern,
                   max.mZ.Da.error = max.mZ.Da.error, 
                   ignoreMascotIonScore = ignoreMascotIonScore,
                   topN = topN,
                   fragmentIonMzRange = fragmentIonMzRange,
                   fragmentIonRange = fragmentIonRange),
        rt.normalized=unlist(x.rt), 
        rt.input=unlist(lapply(data, function(x){x$rt}))))
}

