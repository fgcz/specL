#R
# $HeadURL$
# $Id$
# $Date$


test_genSwathIonLib<-
function(){
# TEST 1
    r.genSwathIonLib.top5 <- genSwathIonLib(peptideStd,
    peptideStd.redundant, topN=5, 
    fragmentIonRange=c(2,100),
    fragmentIonMzRange=c(200,2000),
    fragmentIonFUN=function (b, y) {
      return( cbind(y1_=y) )
      }
    )
    
    checkEqualsNumeric(r.genSwathIonLib.top5@ionlibrary[[40]]@q3, 
    c(800.4497, 604.3285,1016.5222, 503.2805, 929.4925),
    tolerance=0.001)

    checkEqualsNumeric(r.genSwathIonLib.top5@rt.input[1:5], 
                       c(16.83185, 13.13262, 18.54058, 18.36923, 15.30478), 
                       tolerance=0.001)
    
    checkEqualsNumeric(r.genSwathIonLib.top5@rt.normalized[1:5], 
                       c(95.97314,  52.60417, 116.00582, 113.99703,  78.07007), 
                       tolerance=0.001)

    checkEqualsNumeric(length(r.genSwathIonLib.top5@rt.normalized), 
                       length(r.genSwathIonLib.top5@rt.normalized), 
                       tolerance=0.001)
}

test_genSwathIonLib_noiRT_peptides<-
function(){
    # case 1: useless peptides

    #reverse n iRT peptides
    n<-7
    peptide <- unlist(lapply(as.character(iRTpeptides$peptide[1:n]),
                             function(x){paste(rev(strsplit(x,"")[[1]]), collapse='')})) 
    
    rt <- seq(-24.9200, 122.2462, length=length(peptide))
    myiRTpeptides<-cbind(peptide=peptide, rt=rt)

    peptideStd.ionLib <- genSwathIonLib(data=peptideStd, 
                                        data.fit=peptideStd.redundant, 
                                        iRT=myiRTpeptides)

    # 
    res <- lapply (1:length(peptideStd), 
        function(idx){ 
            checkEqualsNumeric(ionlibrary(peptideStd.ionLib)[[idx]]@irt, 
                peptideStd[[idx]]$rt, 
                tolerance=0.05)
            }
            )
}


test_genSwathIonLib()
test_genSwathIonLib_noiRT_peptides()

#TODO(cp):
# check lower boundary of transitions
# genSwathIonLib(peptideStd, peptideStd.redundant, topN=15, fragmentIonRange=c(12, 20))
# nur alle peptide die mind. 12 assigneed fragment ions max 15
