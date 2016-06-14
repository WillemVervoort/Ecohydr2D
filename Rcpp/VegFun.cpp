#include <Rcpp.h>
#include "Soilfun.cpp"
using namespace Rcpp;

// first redefine vegpar
List Veg(Character vtype, List soilpar) {
  // general definitions
  double E_w = 0.01;
  double k = 0.5;
  
  if (vtype == "Grass") {
    double Zr = 50.0;
    double delta = 0.1;
    double psi_ss = -0.09;
    // calculate s_star
    double s_star = pow(psi_ss/soilpar["psi_s_bar"],(-1.0/soilpar["b"]));
    double psi_sw = -3.0;
    
    // calculate s_w
    double s_w = pow(psi_sw/soilpar["psi_s_bar"],(-1.0/soilpar["b"]))
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
      double E_max = 0.33;
    double LAI = 1.5; //Roberts et al. 2000 assuming Kc = 0.8
    double c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
    double fr = 1.0;
    double Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
    double c1 = NULL;
    double fs = NULL;
    Logical DR = FALSE;
    
    // vegmodelling parameters taken out
    
  }
  if (vtype == "TreesDR") {
    double Zr = 100.0;
    double delta = 0.2;
    double psi_ss = -0.12; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
    double s_star = pow(psi_ss/soilpar["psi_s_bar"],(-1.0/soilpar["b"]));
    double psi_sw = -7.0; //MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
    
    // calculate s_w
    double s_w = pow(psi_sw/soilpar["psi_s_bar"],(-1.0/soilpar["b"]))
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
      double E_max = 0.5;
    double LAI = 0.9; //Roberts et al. 2000 assuming Kc = 0.8
    double c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
    double fr = 1.0;
    double Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
    double c1 = 1.0;
    double fs = 0.25;
    Logical DR = TRUE;
    
    // vegmodelling parameters taken out
    
  }
  
  if (vtype == "TreesNoDR") {
    double Zr = 100.0;
    double delta = 0.2;
    double psi_ss = -0.12; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
    double s_star = pow(psi_ss/soilpar["psi_s_bar"],(-1.0/soilpar["b"]));
    double psi_sw = -7.0; //MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
    
    // calculate s_w
    double s_w = pow(psi_sw/soilpar["psi_s_bar"],(-1.0/soilpar["b"]))
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
      double E_max = 0.5;
    double LAI = 0.9; //Roberts et al. 2000 assuming Kc = 0.8
    double c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
    double fr = 1.0;
    double Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
    double c1 = 1.0;
    double fs = 0.25;
    Logical DR = FALSE;
    
    // vegmodelling parameters taken out
    
  }
  
  if (vtype == "Lignum") {
    // suggested values F. van Ogtrop & W. Vervoort
    // No data available?
    
    double Zr = 100.0;
    double delta = 0.2;
    double psi_ss = -0.12; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
    double s_star = pow(psi_ss/soilpar["psi_s_bar"],(-1.0/soilpar["b"]));
    double psi_sw = -5.0;
    
    // calculate s_w
    double s_w = pow(psi_sw/soilpar["psi_s_bar"],(-1.0/soilpar["b"]))
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
      double E_max = 0.33;
    double LAI = 2.03; //Roberts et al. 2000 assuming Kc = 0.8
    double c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
    double fr = 1.0;
    double Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
    double c1 = NULL;
    double fs = NULL;
    Logical DR = TRUE;
    
    // vegmodelling parameters taken out
    
  }
  
  if (vtype == "Bare") {
    //................................................
    // Bare Soil
    // added for 2D project
    // Willem Vervoort 12/06/13
    //.................................................
    
    double Zr = 25.0;
    double delta = 0;
    double psi_ss = -0.01; //Mpa Table 1, F_I & R-I 2004
    // calculate s_star
    double s_star = pow(psi_ss/soilpar["psi_s_bar"],(-1.0/soilpar["b"]));
    double psi_sw = -1.5; //MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
    
    // calculate s_w
    double s_w = pow(psi_sw/soilpar["psi_s_bar"],(-1.0/soilpar["b"]))
      
      // colonisation parameters taken out
      
      // vegetation ET parameters
      double E_max = 0.33;
    double LAI = 1.0; //Roberts et al. 2000 assuming Kc = 0.8
    double c_T = 0.45; // Tarrawarra dataset par in Teuling and Troch 2005
    double fr = 1.0;
    double Ep = 0.65; // cm/day BOM data average Moree
    
    // new groundwater uptake function
    // defines fraction of roots close to groundwater
    double c1 = NULL;
    double fs = NULL;
    Logical DR = FALSE;
    
    // vegmodelling parameters taken out
    
  }
  
  
}
