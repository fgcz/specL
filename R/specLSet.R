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

.specLSetTest<- function(x, list_q1){
    # todo(cp): check also q3 and peptideModSeq
    idx <- findNN(x@q1, list_q1)

    if (abs(list_q1[idx] - x@q1) < 0.01){
            return(FALSE)
    }

    return (TRUE)
}

setMethod(f="merge.specLSet", signature="specLSet", 
    definition=function(object0, object1){
        #todo(cp): could be an argument
        FUN <- .specLSetTest

        object0_q1 <- unlist(lapply(object0@ionlibrary, function(x){return(x@q1)}))
        idx <- order(object0_q1)

        for (i in 1:length(object1@ionlibrary) ){
            x <- object1@ionlibrary[[i]]

            if ( FUN(x, object0_q1[idx]) ){
                object0@ionlibrary <- c(object0@ionlibrary, x)

                object0@rt.input <- c(object0@rt.input, object1@rt.input[i])
                object0@rt.normalized <- c(object0@rt.normalized, object1@rt.normalized[i])
            }

        }
        return (object0)
    }
)

setMethod(f="summary", signature="specLSet", function(object){
    cat("Summary of a \"specLSet\" object.\n\n")

    cat("\nInput:\n")
    cat("\nParameter:\n")

    mapply (function(x, y){cat(paste("\t",x,'=',y,'\n',sep=''))}, 
        names(object@input.parameter), object@input.parameter)
    
#    cat("\nOutput:\n")

    cat("\nNumber of precursor (q1 and peptideModSeq) = ")
    cat(length((unlist(lapply(slot(object,"ionlibrary"), function(x){paste(x@q1, x@peptideModSeq, sep='_')})))))
    cat("\nNumber of unique precursor (q1 and peptideModSeq) = ")
    cat(length(unique(unlist(lapply(slot(object,"ionlibrary"), function(x){paste(x@q1, x@peptideModSeq, sep='_')})))))

    cat("\nFrequency of number of transitions:")
    t<-table(unlist((lapply(slot(object,"ionlibrary"), function(x){length(paste(x@q1, x@q3, x@peptideModSeq))}))))
    print (t)

    cat("\nNumber of annotated precursor = ")
    cat(sum(unlist(lapply(slot(object,"ionlibrary"), function(x){x@proteinInformation != ''}))))

    cat("\nFilename(s)\n")
    files<-((unique(unlist(lapply(slot(object,"ionlibrary"), function(x){x@filename})))))
    for (f in files){
        cat("\t")
        cat (f)
        cat ("\n")
    }
    cat("\nMisc:\n")
    memsize <-  format(object.size(object), units = "b")
    cat("\nMemory usage\t=\t", memsize, "\n")
})

setMethod(f="plot", signature="specLSet", 
          definition=function(x, ...){
              file <- as.factor(unlist( lapply(x@ionlibrary, function(y){ y@filename }) ))
              op<-par(mfrow=c(2,1))
              plot(x@rt.normalized ~ x@rt.input,
                main='specLSet iRT normalization',
                xlab="input retention time (min)",
                ylab="independent retention time",
                col=file)

              hist(x@rt.normalized)
              par(op)
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

