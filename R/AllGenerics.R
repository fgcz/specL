#R
# AllGenerics.R
setGeneric("write.spectronaut",  
    function(x, ...) standardGeneric("write.spectronaut"))

#setGeneric("derive_q3_mass_shift", 
#    function(x, ...) standardGeneric("derive_mass_shift"))

setGeneric("ionlibrary", 
    function(object) standardGeneric("ionlibrary"))

setGeneric("rt.input", 
    function(object) standardGeneric("rt.input"))

setGeneric("rt.normalized", 
    function(object) standardGeneric("rt.normalized"))

setGeneric("merge.specLSet", 
    function(object0, object1, FUN) standardGeneric("merge.specLSet"))
