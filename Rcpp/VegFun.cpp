#include <Rcpp.h>
using namespace Rcpp;

// first redefine vegpar
// [[Rcpp::export]]
List Veg_cpp(std::string vtype, List soilpar) {
  // general definitions
  double E_w = 0.01;
  double k = 0.5;
  double psi_s_bar = soilpar["psi_s_bar"];
  double b = soilpar["b"];
  
  // default is "Grass"
     double Zr = 50.0;
    double delta = 0.1;
    double psi_ss = -0.09;
    // calculate s_star
    double s_star = pow(psi_ss/psi_s_bar,(-1.0/b));
    double psi_sw = -3.0;
    
    // calculate s_w
    double s_w = pow(psi_sw/psi_s_bar,(-1.0/b));
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
      double E_max = 0.33;
    double LAI = 1.5; //Roberts et al. 2000 assuming Kc = 0.8
    double c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
    double fr = 1.0;
    double Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
    double c1 = 0.0;
    double fs = 0.0;
    bool DR = FALSE;
    
    // vegmodelling parameters taken out
    
  
  if (vtype == "TreesDR") {
     Zr = 100.0;
     delta = 0.2;
     psi_ss = -0.12; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
     s_star = pow(psi_ss/psi_s_bar,(-1.0/b));
     psi_sw = -7.0; //MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
    
    // calculate s_w
     s_w = pow(psi_sw/psi_s_bar,(-1.0/b));
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
       E_max = 0.5;
     LAI = 0.9; //Roberts et al. 2000 assuming Kc = 0.8
     c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
     fr = 1.0;
     Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
     c1 = 1.0;
     fs = 0.25;
     DR = TRUE;
    
    // vegmodelling parameters taken out
    
  }
  
  if (vtype == "TreesNoDR") {
     Zr = 100.0;
     delta = 0.2;
     psi_ss = -0.12; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
     s_star = pow(psi_ss/psi_s_bar,(-1.0/b));
     psi_sw = -7.0; //MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
    
    // calculate s_w
     s_w = pow(psi_sw/psi_s_bar,(-1.0/b));
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
       E_max = 0.5;
     LAI = 0.9; //Roberts et al. 2000 assuming Kc = 0.8
     c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
     fr = 1.0;
     Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
     c1 = 1.0;
     fs = 0.25;
     DR = FALSE;
    
    // vegmodelling parameters taken out
    
  }
  
  if (vtype == "Lignum") {
    // suggested values F. van Ogtrop & W. Vervoort
    // No data available?
    
     Zr = 100.0;
     delta = 0.2;
     psi_ss = -0.12; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
     s_star = pow(psi_ss/psi_s_bar,(-1.0/b));
     psi_sw = -5.0;
    
    // calculate s_w
     s_w = pow(psi_sw/psi_s_bar,(-1.0/b));
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
       E_max = 0.33;
     LAI = 2.03; //Roberts et al. 2000 assuming Kc = 0.8
     c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
     fr = 1.0;
     Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
     c1 = 0.0;
     fs = 0.0;
     DR = TRUE;
    
    // vegmodelling parameters taken out
    
  }
  
  if (vtype == "Bare") {
    //................................................
    // Bare Soil
    // added for 2D project
    // Willem Vervoort 12/06/13
    //.................................................
    
     Zr = 25.0;
     delta = 0;
     psi_ss = -0.01; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
     s_star = pow(psi_ss/psi_s_bar,(-1.0/b));
     psi_sw = -1.5; //MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
    
    // calculate s_w
     s_w = pow(psi_sw/psi_s_bar,(-1.0/b));
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
       E_max = 0.33;
     LAI = 1.0; //Roberts et al. 2000 assuming Kc = 0.8
     c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
     fr = 1.0;
     Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
     c1 = 0.0;
     fs = 0.0;
     DR = FALSE;
    
    // vegmodelling parameters taken out
    
  }
 return(Rcpp::List::create(Rcpp::Named("Zr") = Zr,
                            Rcpp::Named("delta") = delta,
                            Rcpp::Named("s_star") = s_star,
                            Rcpp::Named("s_w") = s_w,
                            Rcpp::Named("E_max") = E_max,
                            Rcpp::Named("LAI") = LAI,
                            Rcpp::Named("E_w") = E_w,
                            Rcpp::Named("k") = k,
                            Rcpp::Named("c_T") = c_T,
                            Rcpp::Named("fr") = fr,
                            Rcpp::Named("Ep") = Ep,
                            Rcpp::Named("c1") = c1,
                            Rcpp::Named("fs") = fs,
                            Rcpp::Named("DR") = DR));
  
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// sourceCpp("soilfun.cpp")
// soilpar <- Soil_cpp("S Clay Loam")
// Veg_cpp("TreesDR", soilpar)
// */
