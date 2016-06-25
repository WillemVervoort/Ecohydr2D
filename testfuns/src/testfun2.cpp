#include <Rcpp.h>
#include "testfun1.hpp"
using namespace Rcpp;

// //simple test function
// // [[Rcpp::export]]
// NumericVector timesTwo(NumericVector x) {
//   return x * 2;
// }

// [[Rcpp::export]]
NumericVector squared(NumericVector x) {
  NumericVector x1 = timesTwo(x);
  NumericVector y = pow(x1,2);
  return y;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
 squared(2)
*/
