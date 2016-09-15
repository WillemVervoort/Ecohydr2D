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
  bool funswitch = Rcpp::as<bool>(vegpar_in["DR"]);
  
  // Change Rainfall into infiltration
  double Phi = In;
  double surfoff = 0.0;
  double incr = deltat;
  float intincr = 1/incr;
  //Rcpp::Rcout << "rest" << intincr << std::endl;
  
  // check if soil bucket full, create runoff
  if(Phi > (1 - last_t_soil_sat)){
    Phi = 1 - last_t_soil_sat;
    surfoff = (In - Phi)*n*Zr;
  } else {
    Phi = In;
    surfoff = 0.0;
  }
  //Rcpp::Rcout << "funswitch" << funswitch << std::endl;
  NumericVector loss(5); 
  // Call the loss model
  if (funswitch==0) {
    loss = rho_new_cpp(last_t_soil_sat, 
                       soilpar_in, vegpar_in,
                       Z_in,Zmean_in, Z_prev);
   // Rcpp::Rcout << "loss0 = " << loss[0] << std::endl;
  } else {
    loss = FB_new_cpp(last_t_soil_sat, 
                      soilpar_in, vegpar_in,
                      Z_in, Zmean_in, Z_prev);
  }
  // update the water balance
  // ss.t is temporary ss
  double ss_t = last_t_soil_sat + Phi - loss[0]/(n*Zr)*intincr;
  //Rcpp::Rcout << "ss_t" << ss_t << std::endl;
  
  // build in safety if capillary flux fills soil with dynamic water table
  if (ss_t > 1) {  // very wet conditions
    soil_sat = 1;
    qcap = (loss[3] - (ss_t - 1)*n*Zr/intincr)*intincr;
  } else { // normal conditions
    soil_sat = ss_t;
    qcap = loss[3]*intincr;
  }
//Rcpp::Rcout << "soil_sat" << soil_sat << std::endl;
  
  static_stress = zeta(soil_sat, ss, sw, q);
  Tg = loss[2]*intincr;
  // Rcpp::Rcout << loss[2];
  Ts = loss[1]*intincr;
  Tt = Tg + Ts;
  if (loss[4] > 0.0) {
    Leakage = loss[4]*intincr;
  }  else Leakage = 0.0;
  //increase with difference between originally calculated and new (pot.reduced) capillary flux
  smloss = (loss[0]*intincr  - (qcap - loss[3]*intincr)); 
  GWrech = Leakage-Tg-qcap;
  //Rcpp::Rcout << GWrech << std::endl;
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
NumericMatrix WBEcoHyd(int t, double R, double ET_in,
                       CharacterVector vtype, List soilpar,
                       NumericVector s_init,
                       NumericVector Zmean,
                       NumericVector GWdepths,
                       NumericVector GWdepths_prev,
                       int deltat,
                       int NX, int NY){
  
  // both R and Etp are single vectors of the rainfall and the ETpotential
  // vtype is a vector of vegetation types along the grid
  // soilpar is the specification of the soiltype (currently 1 type across the grid)
  // s_init are the inital soil moisture levels for the day
  // deltat is the increment in the integration
  
  // Zmean is a vector of long term mean gw depths
  // GWdepths is the gw depths for each gridcell at t
  // GWdepths_prev are the previous gw depths for each gridcell (at t-1)
  
  // define variables
  int m = NX*NY;
  double R_in = 0.0;
  
  // storage output 
  //Matrices
  // this first matrix stores each output (10) from Waterbalancefunction 
  // for single gridcell, but for all sub timesteps
  NumericMatrix Storage_Subday(10,deltat);
  // This stores each daily average output (10 in rows) 
  // for each day for all gridcells
  NumericMatrix Storage_Grid(10,m);
  
  // run through grid
  for (int j = 0;j < m;j++) {
    // Call the vegpar function
    std::string veg_in = as<std::string>(vtype[j]);
    //Rcpp::Rcout << veg_in << std::endl; 
    List vegpar = Veg_cpp(veg_in, soilpar);
    vegpar["Ep"] = ET_in;
    // define soil and vegparameters
    double Zr = vegpar["Zr"];
    double n = soilpar["n"];
    // Now run integration loop
    double s_old = s_init(j);
    for (int p = 0;p < deltat;p++) {
      if (p == 0) {
        // needs to be t - 1 as counters in C++ start at 0
        R_in = R/(n*Zr);
      } else R_in = 0.0;
      //Rcpp::Rcout << R_in << std::endl; 
      
      // run the water balance
      if (p > 0 ) {
        s_old = Storage_Subday(0,(p-1)); 
      } else s_old = s_init(j);
      
      //     // Storage_Subday is not defined
      List Water_out = WB_fun_cpp(vegpar, R_in, s_old,
                                  soilpar, Zmean(j), deltat, 
                                  GWdepths_prev(j), GWdepths(j));
      //Write water balance list output to subdaily vector
      int xsize =  Storage_Subday.nrow();
      for (int k = 0; k < xsize; k++) {
        Storage_Subday(k,p) = Rcpp::as<double>(Water_out(k));
      }
//       if (t > 190 & t < 192) {
          // Rcpp::Rcout << p << std::endl;
//       }
      //
      // store s_init for gridcell till next day
      if (p == 11) {s_init(j) = Storage_Subday(0,p);
      //
      }
    } // close p loop
    // recalculate to daily output and write away
    // produce different output, produce daily values for each grid cell
    int kk = Storage_Grid.nrow();
    //Rcpp::Rcout << kk << std::endl;
    for (int k = 0; k < kk; k++) {
      Storage_Grid(k,j) = sum(Storage_Subday(k,_));
      if (k == 0 || k == 8) Storage_Grid(k,j) = mean(Storage_Subday(k,_));
    }
    //Rcpp::Rcout << j << std::endl;
    
  } // close j loop
  // x = 4;
  // y = 5;
  // NumericMatrix out(x,y);
  return Storage_Grid;
}
