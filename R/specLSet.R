#R 
#
# $HeadURL$
# $Id$
# $Date$

# this function is for normalizing the rt on data
# for building the model data.fit  is used

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
    definition=function(object0, object1){
      
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
    cat("\nNumber of unique precursor (q1 and peptideModSeq) = ")
    cat(length(unique(unlist(lapply(slot(object,"ionlibrary"), function(x){paste(x@q1, x@peptideModSeq, sep='_')})))))

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

    cat("\nNumber of precursors in Filename(s)\n")
    t <- table(unlist(lapply(object@ionlibrary, function(x){x@filename})))

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

setMethod(f="plot", signature="specLSet", 
          definition=function(x, ...){
              file <- as.factor(unlist( lapply(x@ionlibrary, function(y){ y@filename }) ))
              peptide <- unlist( lapply(x@ionlibrary, function(y){ y@peptide_sequence }) )
            
              
              plot(x@rt.normalized ~ x@rt.input,
                   main='specLSet iRT normalization',
                   xlab="input retention time (min)",
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
              
              
              frg <- unlist(lapply(x@ionlibrary, function(xx){paste(xx@frg_type, xx@frg_z, "+",sep='')}))
              frg.table <- table(frg)
              cm<-rainbow(length(frg.table))
              hist(x@rt.normalized)
              
              plot(unlist(lapply(x@ionlibrary, function(xx){rep(xx@irt, length(xx@q3))})), 
                   unlist(lapply(x@ionlibrary, function(xx){xx@q3})), 
                   col=cm[as.factor(frg)],
                   xlab='rt.normalized',
                   ylab='fragment ion mass',
                   pch='_',
                   main='in-silico rt-fragment ion map')
              legend('topleft', names(frg.table), pch=22, col=cm, cex=0.75)

              
              
              barplot(frg.table, col=cm,
                   main='ion type / charge state')
              #legend('topleft', 
              #       paste(names(frg.table), frg.table, sep'='),
              #       pch=22, col=cm, cex=0.5)
                     
             
          })
            

setMethod(f="write.spectronaut", signature="specLSet", 
          definition=function(x, file="specL.txt", ...){
    if (file.exists(file)){
         warning("file already exists. not writting!")
         return()
    }

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

