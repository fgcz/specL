#R
# $HeadURL$
# $Id$
# $Date$

test_read.bibliospec <-
function(){
# TEST 1
    S.peptideStd <- read.bibliospec(system.file("extdata",  name='peptideStd.sqlite', package = "specL"))
    S.peptideStd.redundant <- read.bibliospec(system.file("extdata",  name='peptideStd_redundant.sqlite', package = "specL"))

    r.genSwathIonLib.top5 <- genSwathIonLib(S.peptideStd,
        S.peptideStd.redundant, topN=5, 
        fragmentIonRange=c(2,100),
        fragmentIonMzRange=c(200,2000),
        fragmentIonFUN=function (b, y) {
            return( cbind(y1_=y) )
        }
    )
    
    checkEqualsNumeric(r.genSwathIonLib.top5@ionlibrary[[40]]@q3, 
        c(800.4497, 604.3285,1016.5222, 503.2805, 929.4925),
        tolerance = 0.001)

    checkEqualsNumeric(r.genSwathIonLib.top5@rt.input[1:5], 
        c(16.83185, 13.13262, 18.54058, 18.36923, 15.30478), 
        tolerance = 0.001)
    
    checkEqualsNumeric(r.genSwathIonLib.top5@rt.normalized[1:5], 
        c(95.97314, 52.60417, 116.00582, 113.99703,  78.07007), 
        tolerance = 0.001)

    checkEqualsNumeric(length(r.genSwathIonLib.top5@rt.normalized), 
        length(r.genSwathIonLib.top5@rt.normalized), 
        tolerance = 0.001)
}

test_read.bibliospec()

