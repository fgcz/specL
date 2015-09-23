#include <Rcpp.h>

#include <algorithm>
#include <cmath>
//#include <function>

/*

$HeadURL$
$Id$

*/

class MyComparator : public std::binary_function<double, double,bool> {
  public : 
    bool operator()(double a , double b){
    return (a <= b);
  }
};

RcppExport SEXP lower_bound_(SEXP q, SEXP vec)
{

    Rcpp::NumericVector xq(q);
    Rcpp::NumericVector xvec(vec);
    Rcpp::IntegerVector IDX;

    MyComparator comparator;
    try {

      for (int i =0; i < xq.size(); ++i) {
      //for (int x:xq){
	    IDX.push_back(std::distance(xvec.begin(),
					std::lower_bound(xvec.begin(),
							 xvec.end(), 
							 xq[i],
							 //x,
							 comparator
					)
			  )
		);
	}
	return (IDX);
    }
    catch( ...) {
	::Rf_error("c++ exception (unknown reason)");
    }

    // TODO(cp): handle possible exceptions
}
