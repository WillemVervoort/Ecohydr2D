#ifndef __testfun1__
#define __testfun1__


#include <Rcpp.h>
using namespace Rcpp;

//simple test function
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
#endif //__testfuns__