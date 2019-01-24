#R


# $HeadURL$
# $Id$
# $Date$

# compute dynamic swath windows
cdsw <- function(x, n = 20, overlap = 1.0, ...) {
  if (is.psmSet(x)) {
    x <- unlist(lapply(x, function(x) {
      x$pepmass
    }))
  } else if (class(x) == 'specLSet') {
    x <- unlist(lapply(x@ionlibrary, function(xx) {
      xx@q1
    }))
  }
  # x should be numeric
  if (!is.numeric(x)) {
    warning("can not compute quantils. 'x' is not numeric.")
    return (NULL)
  }
  
  q <- quantile(x, seq(0, 1, length = n + 1))
  q[1] <- q[1] - 0.5
  q[length(q)] <- q[length(q)] + 0.5
  q <- round(q)
  idx <- 1:n
  from <- q[idx] - overlap * 0.5
  
  to <- q[idx + 1] + overlap * 0.5
  width <- 0.5 * (to - from)
  mid <- from + width
  h <- hist(x,
            breaks = q, ...)
  
  data.frame(from, to, mid, width, counts = h$counts)
}
