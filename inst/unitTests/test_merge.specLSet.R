#R

# $HeadURL$
# $Id$
# $Date$

test_merge.specLSet<-
function(){
  
  s1 <- genSwathIonLib(peptideStd[seq(1,length(peptideStd), by=2)], peptideStd.redundant)
  s2 <- genSwathIonLib(peptideStd[seq(2,length(peptideStd), by=2)], peptideStd.redundant)
  s <- genSwathIonLib(peptideStd, peptideStd.redundant)
  
  m12 <- merge.specLSet(s1, s2)
  m21 <- merge.specLSet(s2, s1)

  idx.s <- order(unlist(lapply(ionlibrary(s), function(x){x@group_id})))
  idx.m12 <- order(unlist(lapply(ionlibrary(m12), function(x){x@group_id})))
  idx.m21 <- order(unlist(lapply(ionlibrary(m21), function(x){x@group_id})))
  
  
  checkEqualsNumeric(cor(unlist(lapply(ionlibrary(m12)[idx.m12], function(x){x@q3})), 
      unlist(lapply(ionlibrary(m21)[idx.m21], function(x){x@q3}))), 1, tolerance=0.0)
      
      
  checkEqualsNumeric(cor(unlist(lapply(ionlibrary(s)[idx.s], function(x){x@q3})), 
          unlist(lapply(ionlibrary(m21)[idx.m21], function(x){x@q3}))), 1, tolerance=0.0)
}

test_merge.specLSet()

