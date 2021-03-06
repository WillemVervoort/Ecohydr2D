// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
#include "VegFun.hpp"
#include "Fluxfunctions.hpp"
using namespace Numer;
using namespace Rcpp;

// static stress function
// [[Rcpp::export]]
double zeta(double s,double s_star, double s_w,double q) {
  double temp = std::max(pow((s_star-s)/(s_star-s_w),q),0.0);
  double out = std::min(temp,1.0);
  return out;
}


// Rewriting the 2-D water balance implementation across grid cells
// First rewrite the water balance function
// [[Rcpp::export]]
List WB_fun_cpp(List vegpar_in, double In, double last_t_soil_sat, 
                List soilpar_in, 
                double Zmean_in,
                int deltat,
                double Z_prev,
                double Z_in) 
{
  const double n = soilpar_in["n"];
  const double Zr = vegpar_in["Zr"];
  const double ss = vegpar_in["s_star"];
  const double sw = vegpar_in["s_w"];
  const double q = vegpar_in["q"];
  double soil_sat = 0.0;
  double qcap = 0.0;
  double static_stress = 0.0;
  double Tg = 0.0;
  double Ts = 0.0;
  double Tt = 0.0;
  double Leakage = 0.0;
  double smloss = 0.0;
  double GWrech = 0.0;
  const bool funswitch = vegpar_in["DR"];
  
  // Change Rainfall into infiltration
  double Phi = In;
  double surfoff = 0.0;
  double intincr = 1/deltat;
  // check if soil bucket full, create runoff
  if(Phi > (1 - last_t_soil_sat)){
    Phi = 1 - last_t_soil_sat;
    surfoff = (In - Phi)*n*Zr;
  } 
  NumericVector loss(5); 
  // Call the loss model
  if (funswitch==TRUE) {
    loss = FB_new_cpp(last_t_soil_sat, 
                      soilpar_in, vegpar_in,
                      Z_in, Zmean_in, Z_prev);
    //Rcpp::Rcout << loss[2]
  } else {
    loss = rho_new_cpp(last_t_soil_sat, Z_in,
                       soilpar_in, vegpar_in,
                       Zmean_in, Z_prev);
  }
  // update the water balance
  // ss.t is temporary ss
  double ss_t = last_t_soil_sat + Phi - loss[0]/(n*Zr)*intincr;
  
  // build in safety if capillary flux fills soil with dynamic water table
  if (ss_t > 1) {  // very wet conditions
    soil_sat = 1;
    qcap = (loss[3] - (ss_t - 1)*n*Zr/intincr)*intincr;
  } else { // normal conditions
    soil_sat = ss_t;
    qcap = loss[3]*(1/intincr);
  }
  static_stress = zeta(soil_sat, ss, sw, q);
  Tg = loss[2]*intincr;
  // Rcpp::Rcout << loss[2];
  Ts = loss[1]*intincr;
  Tt = Tg + Ts;
  if (loss[4] > 0.0) {
    Leakage = loss[4]*intincr;
  } 
  //increase with difference between originally calculated and new (pot.reduced) capillary flux
  smloss = (loss[0]*intincr  - (qcap - loss[3]*intincr)); 
  GWrech = Leakage-loss[2]*intincr-qcap;
  //Rcpp::Rcout << qcap;
  List out = Rcpp::List::create(Rcpp::Named("s") = soil_sat,
                                Rcpp::Named("qcap") = qcap, 
                                Rcpp::Named("Tgw") = Tg, 
                                Rcpp::Named("Tsoil") = Ts, 
                                Rcpp::Named("Ttotal") = Tt, 
                                Rcpp::Named("Leakage") = Leakage, 
                                Rcpp::Named("GWrech") = GWrech, 
                                Rcpp::Named("smloss") = smloss, 
                                Rcpp::Named("static_stress") = static_stress,
                                Rcpp::Named("surfoff") = surfoff);
  return out;      
}



// this is in the middle part of the overall Bigfun

// [[Rcpp::export]]
NumericMatrix WBEcoHyd(int x, int y, int t, NumericVector R, NumericVector ETp,
                       CharacterVector vtype, List soilpar,
                       NumericVector s_init,
                       double fullday,
                       NumericVector Zmean,
                       NumericVector GWdepths,
                       NumericVector GWdepths_prev,
                       int deltat,
                       NumericVector NX, NumericVector NY){
  
  // both R and Etp are single vectors of the rainfall and the ETpotential
  // vtype is a vector of vegetation types along the grid
  // soilpar is the specification of the soiltype (currently 1 type across the grid)
  // s_init are the inital soil moisture levels for the day
  // deltat is the increment in the integration
  
  // Zmean is a vector of long term mean gw depths
  // GWdepths is the gw depths for each gridcell at t
  // GWdepths_prev are the previous gw depths for each gridcell (at t-1)
  
  // define variables
  int m = NX.size()*NY.size();
  double R_in = 0.0;
  
  // storage output vectors
  // grid vectors
  NumericVector s_g(m), qcap_g(m), Tgw_g(m), Tsoil_g(m), Ttotal_g(m), Leakage_g(m);
  NumericVector GWrech_g(m), smloss_g(m), static_stress_g(m), surfoff_g(m);
  //Matrices
  NumericMatrix Storage_Subday(12,10);
  // this only works as long as NX = 1
  NumericMatrix Storage_Grid(10,m);
  
  // run through grid
  for (int j = 0;j < m;j++) {
    // Call the vegpar function
    std::string veg_in = as<std::string>(vtype[j]);
    List vegpar = Veg_cpp(veg_in, soilpar);
    vegpar["Ep"] = ETp[j];
    // define soil and vegparameters
    double Zr = vegpar["Zr"];
    double n = soilpar["n"];
    // Now run integration loop
    for (int p = 0;p < deltat;p++) {
      if (p == 0) {
        R_in = R[t]/(n*Zr);
      }
      Rcpp::Rcout << "R_in =" << R_in; 
      // run the water balance
      double s_old = s_init[j];
      if (p > 1 ) s_old = Storage_Subday((p-1),1);
      //     // Storage_Subday is not defined
      List Water_out = WB_fun_cpp(vegpar, R_in, s_old,
                                  soilpar, Zmean[j], deltat, GWdepths_prev[j], GWdepths[j]);
      
      //Write water balance list output to subdaily vector
      int xsize =  Storage_Subday.ncol();
      for (int k = 0; k < xsize; k++) {
        Storage_Subday(p,k) = Rcpp::as<double>(Water_out(k));
      }
      //
      // store s_init for gridcell till next day
      if (p == 12) s_init(j) = Storage_Subday(p,1);
    } // close p loop
    // recalculate to daily output and write away
    // produce different output, produce daily values for each grid cell
    int kk = Storage_Grid.ncol();
    for (int k = 0; k < kk; k++) {
      Storage_Grid(j,k) = sum(Storage_Subday(_,k));
      if (k > 0) Storage_Grid(j,k) = mean(Storage_Subday(_,k));
    }
    
  } // close j loop
  // x = 4;
  // y = 5;
  // NumericMatrix out(x,y);
  return Storage_Grid;
}
