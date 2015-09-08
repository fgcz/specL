#R
.annotateProteinIDGrep <- function(x , fasta, digestPattern ){
  sequence = x$peptideSequence
  idx <- grep (sequence,  fasta, fixed = TRUE)
  
  pattern = paste(digestPattern, sequence, sep='')
  selected <- fasta[idx]
  idx2 <- grep(pattern, selected, fixed=FALSE)
  idx<-idx[idx2]
  x$proteinInformation = names(fasta)[idx]
  return(x)
}

annotate.protein_id <- function(data,  
                                file = NULL,
                                fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA"),
                                digestPattern = "(([RK])|(^)|(^M))"
                                ){
  message("start protein annotation ..." )
  time.start <- Sys.time();
  
    if( length(data) > 200 & parallel::detectCores(logical=FALSE) > 1){
      mc.cores = 1
      if(Sys.info()['sysname'] != "Windows"){
        mc.cores = parallel::detectCores(logical=FALSE)
      }
      data <- parallel::mclapply(data, .annotateProteinIDGrep, fasta, digestPattern, mc.cores =  mc.cores) 
    }else{
      data <- lapply(data, .annotateProteinIDGrep, fasta, digestPattern)
    }
  class(data) <- "psmSet"
  time.end <- Sys.time();
  message(paste("time taken: ", difftime(time.end, time.start, units='mins'),  "minutes"))
  return(data)
}
