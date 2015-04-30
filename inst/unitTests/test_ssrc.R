#R
# $HeadURL$
# $Id$
# $Date$

test_ssrc <-
function(){
        

    ssrc <- unlist(lapply(peptideStd, function(x){ssrc(x$peptideSequence)}))
    rt <- unlist(lapply(peptideStd, function(x){x$rt}))

    checkEqualsNumeric(cor(ssrc, rt, method='spearman'), 0.91, tolerance = 0.1)

}


test_ssrc()
