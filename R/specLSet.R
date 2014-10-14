#R 
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/swath.R $
# $Id: swath.R 6657 2014-09-10 12:21:01Z cpanse $
# $Date: 2014-09-10 14:21:01 +0200 (Wed, 10 Sep 2014) $

# this function is for normalizing the rt on data
# for building the model data.fit  is used

#' An S4 class to represent a specLSet result
#'
#' @slot contains single elements or vectors

specLSet <- setClass("specLSet",
  slots=c(ionlibrary="list",
  rt.normalized="numeric",
  rt.input="numeric"),
)

setMethod(f="show", signature="specLSet", function(object){
    cat("An \"specLSet\" object.\n\n")
    cat("\ncontent:\n")
    lapply(c("rt.normalized", "rt.input"), function(x){cat(x, '=', slot(object,x), '\n', fill=TRUE)})
    lapply(slot(object,"ionlibrary"), show)
    cat("\nsize:\n")
    memsize <-  format(object.size(object), units = "b")
    cat("Memory usage:", memsize, "\n")

})

setMethod(f="plot", signature="specLSet", 
          definition=function(x, ...){
              file <- as.factor(unlist( lapply(x@ionlibrary, function(y){ y@filename }) ))
              plot(x@rt.normalized ~ x@rt.input,
                main='specLSet iRT normalization',
                xlab="input retention time (min)",
                ylab="independent retention time",
                col=file)
          })
            

setMethod(f="write.Spectronaut", signature="specLSet", 
          definition=function(x, file="specL.txt", ...){
    if (file.exists(file)){
         warning("file already exists. not writting!")
         return()
    }

    res <- lapply(x@ionlibrary, function(xx){write.Spectronaut(xx, file=file)})
}) 

setMethod(f="ionlibrary",  signature="specLSet", 
    definition=function(object) object@ionlibrary)

setMethod(f="rt.input",  signature="specLSet",
    definition=function(object) object@rt.input)

setMethod(f="rt.normalized",  signature="specLSet", 
    definition=function(object) object@rt.normalized)

