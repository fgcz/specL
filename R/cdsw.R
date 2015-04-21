#R


# $HeadURL$
# $Id$
# $Date$

# compute dynamic swath windows
cdsw <- function(x, n=20, overlap=1.0){
  if (class(x) == "psmSet") {
    x <- unlist(lapply(x, function(x){x$pepmass}))
  }
  
  # x should be numeric
  if (class(x) != "numeric"){
    warning("can not compute quantules. 'x' is not numeric.")
  }
  
  q <- quantile(x, seq(0, 1, length=n + 1))
  q[1] <- q[1] - 0.5
  q[length(q)] <- q[length(q)] + 0.5
  q <- round(q)
  idx <- 1:n
  from <- q[idx] - overlap * 0.5;
  to <- q[idx+1] + overlap * 0.5
  width <- 0.5 * (to - from)
  mid <- from + width 
  h <- hist(x, breaks=q, plot=FALSE)
  res <- cbind(from, to, mid, width, counts=h$counts)
  return(res)
}