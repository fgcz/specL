#R 
# $HeadURL$
# $Id$
# $Date$

# this function is for normalizing the rt on data
# for building the model data.fit  is used

#' An S4 class to represent a specL result
#'
#' @slot contains single elements or vectors
specL <- setClass("specL",
  slots=c(group_id="character",
    peptide_sequence="character",
    proteinInformation="character",
    q1="numeric",
    q1.in_silico="numeric",
    q3="numeric",
    q3.in_silico="numeric", 
    decoy="character",
    # precursor charge
    prec_z="numeric", 
    frg_type="character", 
    frg_nr="numeric", 
    # charge of fragment
    frg_z="numeric", 
    relativeFragmentIntensity="numeric", 
    irt="numeric", 
    peptideModSeq="character", 
    mZ.error="numeric",
    filename="character")
)
#TODO include varMods also in the specL class

setMethod("show", "specL", function(object){
    cat("An \"specL\" object.\n\n")
    cat("\ncontent:\n")
    lapply(slotNames(object), function(x){
      v <- slot(object, x)
      if (x == 'filename'){
        cat (substr(v, nchar(v)-getOption("width")-1, nchar(v)))
      }else{
      cat(x, '=', slot(object,x), '', fill=TRUE)}
      }
      )

    cat("\nsize:\n")
    memsize <-  format(object.size(object), units = "b")
    cat("Memory usage:", memsize, "\n")

})

setMethod(f="plot", signature="specL", 
          definition=function(x, y, ...){
  q3<-slot(x, "q3")
  relativeFragmentIntensity <- slot(x, "relativeFragmentIntensity")
  #n<-nchar(as.character(unique(file)))
  #filename.short <- substr(as.character(unique(x@filename)), n - 25, n)
  
  plot(relativeFragmentIntensity ~ q3, 
     main=slot(x, "group_id"),
     sub=paste(x@proteinInformation),
     ylim=c(1,110),
     col='#AAAAAAAA',
     type='h',
     ...)

  points(q3, relativeFragmentIntensity, pch=22)
  text(q3, relativeFragmentIntensity, 
       paste(slot(x, "frg_type"), slot(x, "frg_nr"), sep="_" ),
       col='red', pos=3)
  
  text(q3, relativeFragmentIntensity, 
       paste("z", slot(x, "frg_z"), sep='='),
       col='blue', pos=1)
  
  
  legend("topleft", 
         c(paste("q1", round(slot(x,"q1"),2)), 
                 paste("irt", round(slot(x,"irt"),2))
                 ))  
})

#TODO(cp):
#setMethod(f="derive_q3_mass_shift", signature="specL", 
#          definition=function(x, 
#            shift=list(AA=c('R', 'K'), deltamass=c(-10.008269, -8.014199)), ...){ 
#
#        # S4 copy constructor?
#        y <- x
#
#        for (i in 1:length(shift$AA)){
#
#            lastAA <- tail(strsplit(x@peptide_sequence, '')[[1]], n=1)
#
#            if (lastAA == shift$AA[i]){
#
#                y@q3[x@frg_type == 'y'] <- x@q3[x@frg_type == 'y'] + shift$deltamass[i]
#
#                break
#
#            }
##        }
#
#    return (y)
#}
#)

setMethod(f="write.spectronaut", signature="specL", 
          definition=function(x, file="specL.txt",protIDSeparator=";", ...){
            
            data=cbind(group_id=x@group_id,
                       peptide_sequence=x@peptide_sequence,
                       q1=x@q1,
                       q3=x@q3,
                       q3.in_silico=x@q3.in_silico,
                       decoy=x@decoy,
                       prec_z=x@prec_z,
                       frg_type=x@frg_type,
                       frg_nr=x@frg_nr,
                       frg_z=x@frg_z,
                       relativeFragmentIntensity=x@relativeFragmentIntensity,
                       irt_or_rt=x@irt,
                       peptideModSeq=x@peptideModSeq,
                       mZ.error=x@mZ.error,
                       proteinInformation = paste(x@proteinInformation,collapse=protIDSeparator),
                       filename=x@filename
            )
            
            if (file.exists(file) == FALSE){
              message(paste("writting specL object (including header) to file '", file,  "' ...", sep=''))

              write.table(data, file=file, row.names = FALSE, 
                          col.names = TRUE, quote = FALSE, sep = "\t", append=FALSE)    

            }else{

              write.table(data, file=file, row.names = FALSE, 
                          col.names = FALSE, quote = FALSE, sep = "\t", append=TRUE)           

            }
}
)
          
