# include <Rcpp.h>
// this does not work,check out:
// http://stackoverflow.com/questions/13995266/using-3rd-party-header-files-with-rcpp
using namespace Rcpp;


// Rewriting the 2-D water balance implementation across grid cells


// [[Rcpp::export]]
List WBEcoHyd(NumericVector R, NumericVector ETp,
              CharacterVector vtype, List soilpar.m,
              NumericVector last_t_soil_sat,
              double fullday, 
              Character model,
              double Zmean,m,
              double deltat.m,
              double fs,
              NumericVector Z.prev,
              NumericVector Z.m,
              NumericVector NX, NumericVector NY) {
  
  // define variables
  int j 
  int p
  int n_int = 12
  int m = NX.size()*NY.size()
  // storage output vectors
  // sub daily vectors
  NumericVector s(n_int), qcap(n_int), Tgw(n_int), Tsoil(n_int), Ttotal(n_int), Leakage(n_int)
  NumericVector GWrech(n_int), smloss(n_int), static_stress(n_int), surfoff(n_int)
  // grid vectors
  NumericVector s_g(m), qcap_g(m), Tgw_g(m), Tsoil_g(m), Ttotal_g(m), Leakage_g(m)
  NumericVector GWrech_g(m), smloss_g(m), static_stress_g(m), surfoff_g(m)
  
  // run through grid
  for (j, j++, m-1) {
    // Call the vegpar function
    List vegpar = Veg(vtype(j), soilpar.m);
    vegpar["fs"] = fs;
    vegpar["Ep"] = ETp[t,2];
    vegpar["q"] = 1.0;
    
    // Define which model will be used
    if (vegpar["DR"]==TRUE) {
      model_in = "FB_new";
    } else {
      model_in = "rho_new_1";
    }
    // Now run integration loop
    for (p, p++, n_int-1) {
      if (p == 1) {
        double R.in = R[t]]/(soilpar["n"]*soilpar["Zr"]);
      } else {
        double R.in = 0.0;
      }
      
      // run the water balance
      
      
    }
    
    
    //How to call a R function using Rcpp
    NumericVector callFunction(NumericVector x, Function f) {
      NumericVector res = f(x);
      return res;
    }
  }