// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// M(t) = E(exp(t * X)) = int exp(t * x) * f(x) dx, f(x) is the p.d.f.
class test_f: public Func
{
private:
  const double c1;
  const double z;
public:
  test_f(double c1_,  double z_) : c1(c1_), z(z_) {}
  
  double operator()(const double& x) const
  {
    return std::exp(-c1*z/100);
  }
};

// equation 8 VvdZ 2009
// [[Rcpp::export]]
double int_test(double c1, double z,double z1,double z2) {
  test_f f(c1,z);
  double err_est;
  int err_code;
  return integrate(f,z1, z2, err_est, err_code)/integrate(f,0,1E12, err_est, err_code);
}

// new deep root functions following Orellana et al. 2012
double s_fun(double x,double fs) {
 return exp(-(fs*x)); 
}


double  RWU(double z1,double Zmean, double fs) {
  return R::dnorm(z1/100,Zmean/100,s_fun(Zmean/100,fs),FALSE)*2; // test values
} 
// [[Rcpp::export]]
double Rc_B(double z1,double z2,double c1,double Zmean,double fs) {
  return RWU(z1,Zmean,fs)*int_test(c1,1000,z1,z2);
}
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically 
// // run after the compilation.
// //
// 
// 
// 
// /*** R
//   Rc_B(0,100,1,50,2.5)
// */
