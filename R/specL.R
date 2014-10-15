#R 
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/swath.R $
# $Id: swath.R 6657 2014-09-10 12:21:01Z cpanse $
# $Date: 2014-09-10 14:21:01 +0200 (Wed, 10 Sep 2014) $

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
    q3="numeric",
    q3.in_silico="numeric", 
    decoy="character", 
    prec_z="numeric", 
    frg_type="character", 
    frg_nr="numeric", 
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
    lapply(slotNames(object), function(x){cat(x, '=', slot(object,x), '', fill=TRUE)})

    cat("\nsize:\n")
    memsize <-  format(object.size(object), units = "b")
    cat("Memory usage:", memsize, "\n")

})

setMethod(f="plot", signature="specL", 
          definition=function(x, y, ...){
  
  q3<-slot(x, "q3")
  relativeFragmentIntensity <- slot(x, "relativeFragmentIntensity")
  
  plot(relativeFragmentIntensity ~ q3, 
     main=slot(x, "group_id"),
     sub=x@proteinInformation,
     ylim=c(0,110),
     col='#AAAAAAAA',
     type='h')

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



setMethod(f="write.Spectronaut", signature="specL", 
          definition=function(x, file="specL.txt", ...){
            
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
                       irt=x@irt,
                       peptideModSeq=x@peptideModSeq,
                       mZ.error=x@mZ.error,
                       proteinInformation=x@proteinInformation
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
          
