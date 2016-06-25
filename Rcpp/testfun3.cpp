// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Numer;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix test(int x, int y) {
  NumericMatrix out(x,y);
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test(5,4)
*/
