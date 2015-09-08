#R 
#
# $HeadURL$
# $Id$
# $Date$
#
#' An S4 class to represent a specLSet result
#'
#' @slot contains single elements or vectors

specLSet <- setClass("specLSet",
                     slots=c(ionlibrary="list",
                             rt.normalized="numeric",
                             rt.input="numeric",
                             input.parameter="list")
)

setMethod(f="show", signature="specLSet", function(object){
  cat("An \"specLSet\" object.\n\n")
  cat("\ncontent:\n")
  lapply(c("rt.normalized", "rt.input"), 
         function(x){cat(x, '=', slot(object,x), '\n', fill=TRUE)})
  lapply(slot(object,"ionlibrary"), show)
  cat("\nsize:\n")
  memsize <-  format(object.size(object), units = "b")
  cat("Memory usage:", memsize, "\n")
})

setMethod(f="merge.specLSet", signature="specLSet", 
          definition=function(x, y, ... ){
            object0<-x
            object1<-y
            #todo(cp): check parameter of both objects; if one of the parameter differs 
            # raise an exception
            
            object0_group_id <- unlist(lapply(object0@ionlibrary, function(x){return(x@group_id)}))
            
            for (i in 1:length(object1@ionlibrary) ){
              x <- object1@ionlibrary[[i]]
              
              #if ( FUN(x, object0_q1[idx]) ){\
              if (x@group_id %in% object0_group_id){
                
              }else{
                object0@ionlibrary <- c(object0@ionlibrary, x)
                
                object0@rt.input <- c(object0@rt.input, object1@rt.input[i])
                object0@rt.normalized <- c(object0@rt.normalized, object1@rt.normalized[i])
              }
              
            }
            return (object0)
          }
)

setMethod(f="summary", signature="specLSet", function(object){
  cat("Summary of a \"specLSet\" object.\n")
  #
  #    cat("\nInput:\n")
  cat("\nParameter:\n")
  
  mapply (function(x, y){cat(paste("\t",x,'=',y,'\n',sep=''))}, 
          names(object@input.parameter), object@input.parameter)
  
  #    cat("\nOutput:\n")
  
  cat("\nNumber of precursor (q1 and peptideModSeq) = ")
  cat(length((unlist(lapply(slot(object,"ionlibrary"), function(x){paste(x@q1, x@peptideModSeq, sep='_')})))))
  cat("\nNumber of unique precursor\n(q1.in-silico and peptideModSeq) = ")
  cat(length(unique(unlist(lapply(slot(object,"ionlibrary"), function(x){paste(x@q1.in_silico, x@peptideModSeq, sep='_')})))))
  
  cat("\nNumber of iRT peptide(s) = ")
  cat(sum(unlist(lapply(slot(object,"ionlibrary"), function(x){
    if (x@peptide_sequence %in% iRTpeptides$peptide){1}else{0}}))))
  
  cat("\nNumber of transitions frequency:\n")
  t<-table(unlist((lapply(slot(object,"ionlibrary"), function(x){length(paste(x@q1, x@q3, x@peptideModSeq))}))))
  
  n <- names(t)
  w <- getOption("width") / 2
  for (i in 1:length(t)){
    cat('\t')
    cat(substr(n[i], nchar(n[i])-w, nchar(n[i])))
    cat('\t')
    cat(t[i])
    cat('\n')
  }
  
  
  cat("\nNumber of annotated precursor = ")
  cat(sum(unlist(lapply(slot(object,"ionlibrary"), function(x){x@proteinInformation != ''}))))
  
  
  t <- table(unlist(lapply(object@ionlibrary, function(x){x@filename})))
  cat("\nNumber of file(s)\n")
  cat('\t')
  cat(length(t))
  cat('\n')
  cat("\nNumber of precursors in Filename(s)\n")
  n <- names(t)
  w <- getOption("width") / 2
  for (i in 1:length(t)){
    cat('\t')
    cat(substr(n[i], nchar(n[i])-w, nchar(n[i])))
    cat('\t')
    cat(t[i])
    cat('\n')
  }    
  
  cat("\nMisc:\n")
  memsize <-  format(object.size(object), units = "b")
  cat("\nMemory usage\t=\t", memsize, "\n")
})



.retentiontimePlot <- function(x, file, peptide, ...){
  plot(x@rt.normalized ~ x@rt.input,
       main='specLSet iRT normalization',
       xlab="input retention time",
       ylab="independent retention time",
       col=file, ...)
  
  idx <- which(peptide %in% iRTpeptides$peptide)
  # mark the iRT peptides 
  points(x@rt.input[idx], 
         x@rt.normalized[idx], 
         col=file[idx], 
         lwd=4, pch="x", cex=1.5)
  
  legend('topleft', 'iRT peptides', pch='x')
  
  n<-nchar(as.character(unique(file)))
  legend('bottomright', 
         substr(as.character(unique(file)), n - 25, n),
         pch=22, 
         col=unique(file),
         title='input file names',
         cex=1.0)
}

.ionChargeState <- function(x, frg, frgTable){
  cm<-rainbow(length(frgTable))
  hist(x@rt.normalized)
  
  barplot(frgTable, col=cm,
          main='ion type / charge state')
}


setMethod(f="plot", signature="specLSet", 
          definition=function(x, iRTpeptides=specL::iRTpeptides, art=FALSE, ...){
            
            file <- as.factor(unlist( lapply(x@ionlibrary, function(y){ y@filename }) ))
            peptide <- unlist( lapply(x@ionlibrary, function(y){ y@peptide_sequence }) )
            
            .retentiontimePlot(x, file, peptide, ...)
            
            frg <- unlist(lapply(x@ionlibrary, function(xx){paste(xx@frg_type, xx@frg_z, "+",sep='')}))
            frgTable <- table(frg)
            
            .ionChargeState(x,frg, frgTable)
            
            
            
            #legend('topleft', 
            #       paste(names(frgTable), frgTable, sep'='),
            #       pch=22, col=cm, cex=0.5)
            
            # todo(cp): make a 2nd use case if to many data items are avalable
            if (art==TRUE && require(plotrix)){
              colorMapAlpha<-rainbow(length(frgTable), alpha=0.25)
              plot(unlist(lapply(x@ionlibrary, function(xx){rep(xx@q1, length(xx@q3))})), 
                   unlist(lapply(x@ionlibrary, function(xx){xx@q3})), 
                   col=colorMapAlpha[as.factor(frg)],
                   xlab='q1',
                   ylab='q3',
                   main='Relative Fragment Intensity ~ q1 * q3',
                   sub='The areas represent the fragment ion intensities.',
                   type='n',
                   ...)
              
              draw.circle(unlist(lapply(x@ionlibrary, function(xx){rep(xx@q1, length(xx@q3))})),
                          unlist(lapply(x@ionlibrary, function(xx){xx@q3})),
                          unlist(lapply(x@ionlibrary, function(xx){sqrt(xx@relativeFragmentIntensity)})),
                          col=colorMapAlpha[as.factor(frg)], border=colorMapAlpha[as.factor(frg)])
              legend('topleft', names(frgTable), pch=22, col=colorMapAlpha, cex=0.75)
              
              plot(unlist(lapply(x@ionlibrary, function(xx){rep(xx@irt, length(xx@q3))})), 
                   unlist(lapply(x@ionlibrary, function(xx){xx@q3})), 
                   col=colorMapAlpha[as.factor(frg)],
                   xlab='rt.normalized',
                   ylab='fragment ion mass',
                   type='n',
                   main='Relative Fragment Intensity ~ rt * fragment ion map',
                   sub='The areas represent the fragment ion intensities.',
                   ...)
              draw.circle(unlist(lapply(x@ionlibrary, function(xx){rep(xx@irt, length(xx@q3))})),
                          unlist(lapply(x@ionlibrary, function(xx){xx@q3})),
                          unlist(lapply(x@ionlibrary, function(xx){sqrt(xx@relativeFragmentIntensity)}))/10,
                          col=colorMapAlpha[as.factor(frg)], border=colorMapAlpha[as.factor(frg)])
              legend('topleft', names(frgTable), pch=22, col=colorMapAlpha, cex=0.75)
            }
          })


setMethod(f="write.spectronaut", signature="specLSet", 
          definition=function(x, file="specL.txt", ...){
            res <- lapply(x@ionlibrary, function(xx){write.spectronaut(xx, file=file)})
          }) 

#setMethod(f="derive_q3_mass_shift", signature="specLSet", 
#          definition=function(x, 
#            shift=list(AA=c('R', 'K'), deltamass=c(-10.008269, -8.014199)), ...){
#
#    y <- specLSet(lapply(x@ionlibrary, function(xx){derive_q3_mass_shift(xx), shift=shift)}), 
#        rt.normalized=x@rt.normalized,
#        rt.input=x@rt.input)
#    return (y)
#
#}) 

setMethod(f="ionlibrary",  signature="specLSet", 
          definition=function(object) object@ionlibrary)

setMethod(f="rt.input",  signature="specLSet",
          definition=function(object) object@rt.input)

setMethod(f="rt.normalized",  signature="specLSet", 
          definition=function(object) object@rt.normalized)


setMethod(f="generate.consensus",  signature="specLSet", 
          definition=function(object, ...){
            
            id<-unlist(lapply(slot(object,"ionlibrary"), function(x){paste(x@q1.in_silico, x@peptideModSeq, sep='_')}))
            groups <-tapply(1:length(id), id, function(x){x}, simplify=TRUE)
            
            if (sum(unlist(lapply(groups, function(x){length(x)>2}))) == 0){
              warning("no specs for grouping.")
              return (object)
            }
            
            cat("\nCombining specs having the same group_id:\n")
            lapply((which(unlist(lapply(groups, function(x){length(x)>2})))), function(xx){
              cat(paste(groups[[xx]],collapse=', ')[1]);
              cat(paste(" ->",xx,'\n'))})
            cat("\n")
            
            
            output<-lapply(1:length(groups), function(i){
              
              frg_type<-unlist(lapply(groups[[i]], function(x){object@ionlibrary[[x]]@frg_type}))
              fgr_nr<-unlist(lapply(groups[[i]], function(x){object@ionlibrary[[x]]@frg_nr}))
              fgr_z<-unlist(lapply(groups[[i]], function(x){object@ionlibrary[[x]]@frg_z}))
              
              q3<-unlist(lapply(groups[[i]], function(x){object@ionlibrary[[x]]@q3}))
              irt<-unlist(lapply(groups[[i]], function(x){object@ionlibrary[[x]]@irt}))
              relativeFragmentIntensity<-unlist(lapply(groups[[i]], function(x){l<-object@ionlibrary[[x]]; l@relativeFragmentIntensity}))
              
              by.group<-list(frg_type=frg_type,fgr_nr=fgr_nr,fgr_z=fgr_z)
              
              q3 <- aggregate(q3, by=by.group, FUN=mean)
              
              
              relativeFragmentIntensity <- aggregate(relativeFragmentIntensity, by=by.group, FUN=mean)
              
              filename <- unique(sort(unlist(lapply(groups[[i]], function(x){object@ionlibrary[[x]]@filename}))))
              
              if (length(filename)>1){
                filename<-"*consensus*"
              }
              
              res <- specL(group_id=object@ionlibrary[[groups[[i]][1]]]@group_id, 
                           peptide_sequence=object@ionlibrary[[groups[[i]][1]]]@peptide_sequence, 
                           proteinInformation=object@ionlibrary[[groups[[i]][1]]]@proteinInformation,
                           q1=object@ionlibrary[[groups[[i]][1]]]@q1, 
                           q1.in_silico = object@ionlibrary[[groups[[i]][1]]]@q1.in_silico,
                           q3=q3$x, 
                           q3.in_silico=object@ionlibrary[[groups[[i]][1]]]@q3.in_silico, 
                           decoy=object@ionlibrary[[groups[[i]][1]]]@decoy,
                           prec_z=object@ionlibrary[[groups[[i]][1]]]@prec_z, 
                           frg_type=q3$frg_type, 
                           frg_nr=q3[,2], 
                           frg_z=q3[,3], 
                           relativeFragmentIntensity=relativeFragmentIntensity$x, 
                           irt=mean(irt), 
                           peptideModSeq=object@ionlibrary[[groups[[i]][1]]]@peptideModSeq, 
                           mZ.error=rep(-1, nrow(q3)),
                           filename=filename
              )
            })
            
            
            return(specLSet(ionlibrary=output, 
                            input.parameter=list(comment='consensus specLSet object'),
                            rt.normalized=unlist(lapply(output, function(x){x@irt})), 
                            rt.input=unlist(lapply(output, function(x){x@irt}))))
          })


setGeneric("getProteinPeptideTable", 
           function(object,...) standardGeneric("getProteinPeptideTable"))

setMethod(f="getProteinPeptideTable",  signature="specLSet",
          definition=function(object, ...){
            specLibrary = object
            extractPrecursorCharge = function(tmp){
              peps <- c("peptideSequence"=tmp@peptide_sequence,
                        "peptideModSequence" = tmp@peptideModSeq,
                        "z"=tmp@prec_z,
                        "m"=tmp@q1)
              prots <- tmp@proteinInformation
              tmp = NULL
              if(length(prots) > 0){
                for(i in 1:length(prots)){
                  tmp = rbind(tmp, peps)
                }
                return(cbind(prots,tmp))
              }else{
                return(cbind("NA",t(peps)))
              }
              
            }
            res <- lapply(specLibrary@ionlibrary,extractPrecursorCharge)
            res2 <- NULL
            for(i in 1:length(res)){
              res2 <- rbind(res2, res[[i]] )
            }
            return(res2)
          })
