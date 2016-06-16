// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;


class test_f: public Func
{
private:
  const double c1;
public:
  test_f(double c1_) : c1(c1_) {}
  
  double operator()(const double& z) const
  {
    return std::exp(-c1*(z/100));
  }
};

// equation 8 VvdZ 2009
// [[Rcpp::export]]
double int_test(double c1,double z1,double z2) {
  test_f f(c1);
  double err_est;
  int err_code;
  double res = integrate(f,z2,z1, err_est, err_code)/integrate(f,0,10000, err_est, err_code);
  return res;
}



// new deep root functions following Orellana et al. 2012
double s_fun(double x,double fs) {
 return exp(-(fs*x)); 
}

// [[Rcpp::export]]
double  RWU(double z1,double Zmean, double fs) {
  return 2*R::dnorm(z1/100,Zmean/100,s_fun(Zmean/100,fs),FALSE); // test values
}

// [[Rcpp::export]]
double Rc_B_cpp(double z1,double z2,double c1,double Zmean,double fs) {
  return RWU(z1,Zmean,fs)*int_test(c1,z1,z2);
}
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically 
// // run after the compilation.
// //
// 
// 
// 
/*** R
    Rc_B_cpp(0,100,3,100,1.5)
*/
