#R

annotate.protein_id <- function(data,  
    file = NULL,
    fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA"),
    digestPattern = "(([RK])|(^)|(^M))"){
#   ncores = parallel::detectCores()){

    .annotateProteinIDGrep <- function(x){
        idx <- grep (pattern = paste(digestPattern, x$peptideSequence, sep=''),  x = fasta, fixed = FALSE) 
        x$proteinInformation=paste(fasta.names[idx], collapse=";")
        x
    }

    fasta.names <- names(fasta)

    message("start protein annotation ..." )
    time.start <- Sys.time(); 

    if (require(BiocParallel)){
        message("using BiocParallel::bplapply( ..." )
        data <- BiocParallel::bplapply(data, .annotateProteinIDGrep) 
    }else{
        data <- lapply(data, .annotateProteinIDGrep)
    }


    time.end <- Sys.time();
    message(paste("time taken: ", difftime(time.end, time.start, units='mins'),  "minutes"))

    return(data)
}

.annotate.protein.mclapply <- function(data,  
    file = NULL,
    fasta = read.fasta(file = file, as.string = TRUE, seqtype="AA"),
    digestPattern = "(([RK])|(^)|(^M))",
    ncores = parallel::detectCores()){

    .annotateProteinIDGrep <- function(x){
        idx <- grep (pattern = paste(digestPattern, x$peptideSequence, sep=''),  x = fasta, fixed = FALSE) 
        x$proteinInformation=paste(fasta.names[idx], collapse=";")
        x
    }

    fasta.names <- names(fasta)

    message("start protein annotation ..." )
    start.time <- Sys.time(); 

    if ( length(data) > ncores & ncores > 1 ) {

        message(paste ("detected",
            parallel::detectCores(all.tests = FALSE, logical = FALSE), 
            "compute cores. Running", ncores, "jobs at once ..."))

        data <- parallel::mclapply(data, 
            .annotateProteinIDGrep, 
            mc.cores=ncores, mc.preschedule=TRUE)

    } else {

        data <- lapply(data, .annotateProteinIDGrep)

    }

    end.time <- Sys.time();
    message(paste("time taken: ", time.taken <- end.time - start.time, "minutes"))

    return(data)
}
