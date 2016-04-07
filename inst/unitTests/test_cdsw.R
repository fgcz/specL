#R
# $HeadURL$
# $Id$
# $Date$

test_cdsw <-
function(){
        

    peptideStd.specLSet <- genSwathIonLib(data=peptideStd, data.fit=peptideStd.redundant)

    checkEqualsNumeric(sum(cdsw(peptideStd)$to == cdsw(peptideStd.specLSet)$to) == 20, tolerance = 0.0)


}


test_cdsw()
