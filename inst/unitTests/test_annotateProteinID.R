#R


test_annotate.protein_id<-
function(){

    irtFASTAseq <- paste(">zz|ZZ_FGCZCont0260|",
    "iRT_Protein_with_AAAAK_spacers concatenated Biognosys\n",
    "LGGNEQVTRAAAAKGAGSSEPVTGLDAKAAAAKVEATFGVDESNAKAAAAKYILAGVENS",
    "KAAAAKTPVISGGPYEYRAAAAKTPVITGAPYEYRAAAAKDGLDAASYYAPVRAAAAKAD",
    "VTPADFSEWSKAAAAKGTFIIDPGGVIRAAAAKGTFIIDPAAVIRAAAAKLFLQFGAQGS",
    "PFLK\n")

    Tfile <- file();  cat(irtFASTAseq, file = Tfile);
    fasta.irtFASTAseq <-read.fasta(Tfile, as.string=TRUE, seqtype="AA")
    close(Tfile)

    peptideStd <- specL::annotate.protein_id(peptideStd,
        fasta=fasta.irtFASTAseq)

    (idx<-which(unlist(lapply(peptideStd, 
        function(x){nchar(x$proteinInformation)>0}))))


    checkIdentical(unlist(lapply(idx, function(x){peptideStd[[x]]$proteinInformation})), rep("zz|ZZ_FGCZCont0260|",6))

    (idx<-which(unlist(lapply(peptideStd, 
        function(x){x$proteinInformation == "" }))))

    checkEqualsNumeric(sum(idx),  9001, tolerance=0.001)

    checkIdentical(unlist(lapply(idx, function(x){peptideStd[[x]]$proteinInformation})), rep("", 131))
        

}

test_annotate.protein_id()
