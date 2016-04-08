#R
# $HeadURL$
# $Id$
# $Date$

test_cdsw <-
function(){
        

    peptideStd.specLSet <- genSwathIonLib(data=peptideStd, data.fit=peptideStd.redundant)

    checkEqualsNumeric(sum(cdsw(peptideStd, plot = FALSE)$to == cdsw(peptideStd.specLSet, plot = FALSE)$to), 20, tolerance = 0.0)
    checkEqualsNumeric(sum(cdsw(peptideStd, plot = FALSE)$count == cdsw(peptideStd.specLSet, plot = FALSE)$count), 20, tolerance = 0.0)


}


test_cdsw()
