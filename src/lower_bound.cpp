#include <Rcpp.h>

#include <algorithm>
#include <cmath>

/*

$HeadURL$
$Id$

*/

RcppExport SEXP lower_bound_ (SEXP q, SEXP vec) {
  
  Rcpp::NumericVector xq(q);
  Rcpp::NumericVector xvec(vec);
  Rcpp::IntegerVector IDX; 
  
  try{
  
  for(int x : xq){
    IDX.push_back(std::distance(xvec.begin(), 
      std::lower_bound(xvec.begin(), xvec.end(), x, 
        [](const double& a, const double& b){return(a <= b);})
        )
      );
  }
  return(IDX);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  
  // TODO(cp): handle possible exceptions
}

