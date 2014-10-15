#R


test_genSwathIonLib<-
function(){
# TEST 1
    r.genSwathIonLib.top5 <- genSwathIonLib(peptideStd,
    peptideStd.redundant, topN=5,
    fragmentIonFUN=function (b, y) {
      return( cbind(y1_=y) )
      }
    )
    
    checkEqualsNumeric(r.genSwathIonLib.top5@ionlibrary[[40]]@q3, 
    c(800.4497, 604.3285,1016.5222, 503.2805, 929.4925),
    tolerance=0.001)

    checkEqualsNumeric(r.genSwathIonLib.top5@rt.input[1:5], c(16.83185, 13.13262, 18.54058, 18.36923, 15.30478), tolerance=0.001)
    checkEqualsNumeric(r.genSwathIonLib.top5@rt.normalized[1:5], c(95.97314,  52.60417, 116.00582, 113.99703,  78.07007), tolerance=0.001)

    checkEqualsNumeric(length(r.genSwathIonLib.top5@rt.normalized), length(r.genSwathIonLib.top5@rt.normalized), tolerance=0.001)
}


test_genSwathIonLib()
