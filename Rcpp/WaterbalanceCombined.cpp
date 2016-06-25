// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Numer;
using namespace Rcpp;


//
// Soil data as a function
// Ecohydrology WIMEK project
// Willem Vervoort September 2007
// To simplify soil input
// ########################################
// [[Rcpp::export]]
List  Soil_cpp(std::string stype) {
  double psi_sh = -10.0;
  
  // default = Medium Light Clay
  double n = 0.418; // porosity
  // more soil variables for evaporation & losses
  double K_s = 3.51; // cm/day
  double b = 13.48; // neurotheta LMC
  //  double nvg = 1.089;
  // double  avg  = 0.0591;
  double s_fc = 0.364/n; // Field capacity
  double psi_s_bar = -1.5E-3; // This is the bubbling pressure
  double hb = psi_s_bar*(-1e4);
  double spec_y = 0.054; //Johnson 1967 says 0.05, specific yield. 
  
  if (stype == "L Med Clay Stony") {
    // Medium Light Clay
    n = 0.318; // porosity
    // more soil variables for evaporation & losses
    K_s = 3.51; // cm/day
    b = 13.48; // neurotheta LMC
    // nvg = 1.089;
    // avg = 0.0591;
    s_fc = 0.264/n; // Field capacity
    psi_s_bar = -1.5E-3; // This is the bubbling pressure
    hb = psi_s_bar*(-1e4);
    spec_y = 0.054; //Johnson 1967 says 0.05, specific yield. 
  }
  
  if (stype == "S Clay Loam") {
    // Sandy Clay Loam
    n = 0.367; // porosity
    // more soil variables for evaporation & losses
    K_s = 52.08; // cm/day
    b = 6.4069; // neurotheta sandy clay loam
    // avg = 0.0521;
    // nvg = 1.237;
    s_fc = 0.2677/n; // Field capacity
    
    psi_s_bar = -1.2E-3;
    hb = psi_s_bar*(-1e4);
    spec_y = 0.07;  //difference Fc and por Johnson 1967 says 0.07 
  }
  
  if (stype == "Loamy Sand") {
    // Loamy sand
    n = 0.37; // porosity
    // more soil variables for evaporation & losses
    K_s = 175.3; // cm/day
    b = 4.5206;
    // avg = 0.0641;
    // nvg = 1.344;
    s_fc = 0.2098/n; // Field capacity
    
    psi_s_bar = -0.66E-3; // This is the bubbling pressure
    spec_y = 0.17;  // changed to 0.1 to increase rise in gw, not difference por and fc
  }
  
  if (stype == "H Clay") {
    // Medium Heavy Clay
    n = 0.4473; // porosity
    // more soil variables for evaporation & losses
    K_s = 2.82; // cm/day
    b = 16.1501; // neurotheta Medium heavy clay
    // avg = 0.0613;
    // nvg = 1.086;
    s_fc = 0.3936/n; // Field capacity
    
    psi_s_bar = -1.4e-3;
    hb = psi_s_bar*-1.0e4;
    spec_y = 0.05;  // difference por and fc
  }
  
  if (stype == "M Clay") {
    // Medium Clay
    n =0.4391; // porosity
    // more soil variables for evaporation & losses
    K_s = 6.04; // cm/day
    b =  13.5127;
    // avg = 0.0507;
    //doouble nvg = 1.088;
    s_fc = 0.3818/n; // Field capacity
    
    psi_s_bar = -1.75E-3; // This is the bubbling pressure
    hb = psi_s_bar*-1.0e4;
    spec_y = 0.05; // difference por and fc
  }
  
  if (stype == "C Sand") {
    // Coarse Sand
    n = 0.368; // porosity
    // more soil variables for evaporation & losses
    K_s = 182.68; // cm/day
    b = 4.1152;
    // avg = 0.0712;
    // nvg = 1.392;
    s_fc = 0.1895/n; // Field capacity
    
    psi_s_bar = -0.61E-3; // This is the bubbling pressure
    hb = psi_s_bar*-1.0e4;
    spec_y = 0.27;  // difference por and fc
  }
  
  
  // Other derived parameters
  double s_h = pow(psi_sh/psi_s_bar,-1/b); // soil moisture at hygroscopic point
  double beta = 2*b+4; //page 714 Laio et al., 2001a
  
  // Define parameters for Eagleson function
  // Beta parameters
  double  beta1 = 2+3/b;
  // alpha parameter
  double a1 = 1+(3/2)/(beta1-1);
  
  return(Rcpp::List::create(Rcpp::Named("n") = n,
                            Rcpp::Named("K_s") = K_s,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("hb") = hb,
                            Rcpp::Named("psi_s_bar") = psi_s_bar,
                            Rcpp::Named("s_fc") = s_fc,
                            Rcpp::Named("s_h") = s_h,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("beta1") = beta1,
                            Rcpp::Named("a1") = a1,
                            Rcpp::Named("spec_y") = spec_y));
}


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
  double q = 1.0;
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
    q = 1.0;
    
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
    q = 1.0;
    
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
    q = 1.0;
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
    q = 1.0;
    
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
                            Rcpp::Named("DR") = DR,
                            Rcpp::Named("q") = q));
  
  
}

class test_f: public Numer::Func
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
  double res = Numer::integrate(f,z2,z1, err_est, err_code)/Numer::integrate(f,0,10000, err_est, err_code);
  return res;
}



// new deep root functions following Orellana et al. 2012
double s_fun_cpp(double x,double fs) {
  return exp(-(fs*x)); 
}

// [[Rcpp::export]]
double  RWU_cpp(double z1,double Zmean, double fs) {
  return 2*R::dnorm(z1/100,Zmean/100,s_fun_cpp(Zmean/100,fs),FALSE); // test values
}

// [[Rcpp::export]]
double Rc_B_cpp(double z1,double z2,double c1,double Zmean,double fs) {
  return RWU_cpp(z1,Zmean,fs)*int_test(c1,z1,z2);
}


// E_Teuling
// [[Rcpp::export]]
double E_Teuling_cpp(double s, List vegpar) {
  const double sw = vegpar["s_w"]; 
  const double ss = vegpar["s_star"];
  const double fr = vegpar["fr"]; 
  const double c_T = vegpar["c_T"]; 
  const double LAI = vegpar["LAI"]; 
  const double Ep = vegpar["Ep"]; 
  // calculations
  double beta_T = 0.0;
  if(s > sw) {
    if ((s < ss) & (s > sw)) {
      beta_T = (s - sw)/(ss - sw);
    } else {
      beta_T = 1.0;
    }
  }
  double E = fr*beta_T*(1 - exp(-c_T*LAI))*Ep;
  return(E);
}	

// G function
// [[Rcpp::export]]
double G_cpp(double b, double hb, double Z) {
  double b1 = 2.0 + 3.0/b;
  double a1 = 1.0 + (3.0/2.0)/(b1 - 1.0);
  double H1 = a1*pow(hb/Z,b1);
  return(H1);
}

// function to generate m values
// [[Rcpp::export]]
NumericVector m_fun(List vegpar_i, List soilpar_i, double Z_in, double Gin) {
  const double Zr = vegpar_i["Zr"]; 
  const double b = soilpar_i["b"]; 
  const double hb = soilpar_i["hb"]; 
  const double Ks = soilpar_i["K_s"]; 
  const double n = soilpar_i["n"]; 
  const double beta = soilpar_i["beta"]; 
  const double ss = vegpar_i["s_star"];
  
  double s_lim1 = pow((Z_in - Zr)/hb,(-1/b));
  double m = Ks/(n*Zr*(std::exp(beta*(1-s_lim1))-1));
  double m1 = Ks*Gin/(n*Zr*(1-std::exp(beta*(ss-s_lim1))));
  double m2 = Ks*Gin/(n*Zr);
  return(NumericVector::create(m,m1,m2));
}

// rho_new_1
// [[Rcpp::export]]
NumericVector rho_new_cpp(double s, double ZZ, List soilpar, List vegpar,
                          double Z_mean, double Z_prev) {
  // define variables
  double Zr = vegpar["Zr"];
  double ss = vegpar["s_star"];
  double sw = vegpar["s_w"]; 
  double n = soilpar["n"];
  double b = soilpar["b"];
  double hb = soilpar["hb"];
  double beta = soilpar["beta"]; 
  double s_lim = pow(((ZZ - Zr)/hb),(-1/b));
  // apply G function
  double G1 = G_cpp(b,hb,(ZZ-Zr));
  // calculate parameters
  NumericVector m_values = m_fun(vegpar, soilpar, ZZ, G1);
  double E_max = E_Teuling_cpp(s=1,vegpar);
  double eta = E_max/(n*Zr);
  double Q = 0.0;
  double E = 0.0;
  
  // Now calculate the soil moisture
  if (s > s_lim) {
    Q = m_values[0]*(exp(beta*(s - s_lim))-1);
    E = eta;
  }
  if (m_values[1] < eta) {
    // define s_cr: s-critical
    double s_cr = m_values[3]/eta*(ss - sw) + sw;
    
    if (s < s_cr) {
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G_cpp(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      E = 0.0;  
      Q = -(m_values[2]*((s_cr-sw)/(ss-sw))-m_values_prev[2]*((s-sw)/(ss-sw)));
    }
    if (s > s_cr && s <= ss)
    {
      E = eta*((s-s_cr)/(ss-s_cr));
      Q = -m_values[2]*((s-s_cr)/(ss - s_cr));
    } else {
      if (s > ss && s <= s_lim) {
        Q = -m_values[1]*(1-exp(beta*(s-s_lim)));
        E = eta;
      }
    }
  } else {
    if (s > ss && s <= s_lim) {
      E = eta;
      Q = -m_values[1]*exp(beta*(s-s_lim));
    } else {    // again, if GW table rises: E=0 and qcap =difference qcap and old transpiration
      double s_cr = m_values[2]/eta*(ss - sw) + sw;
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G_cpp(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      E = 0.0;  
      Q = -(m_values[2]*((s_cr-sw)/(ss-sw))-m_values_prev[2]*((s-sw)/(ss-sw)));
    }
  }
  double Q_out = Q*(n*Zr);   
  double E_out = E*(n*Zr);   
  double L_out = Q_out + E_out;
  double qcap = 0.0;
  if(Q_out < 0) qcap = abs(Q);
  return(NumericVector::create(L_out,E_out,0.0,qcap,Q_out));
}

// Eagleson function
// [[Rcpp::export]]
NumericVector Fun_E_cpp(double s,double Z,List soilpar,List vegpar) {
  double Zr = vegpar["Zr"];
  double b = soilpar["b"];
  double Ks = soilpar["K_s"];
  double hb = soilpar["hb"];
  double H1 = G_cpp(b,hb,(Z-Zr));
  double Q = Ks*(pow(s,(2*b+3))- H1); // bottom flux
  double  E = E_Teuling_cpp(s,vegpar);
  double qcap = 0.0;
  if (Q < 0) qcap = abs(Q); // capillary flow component
  double L = Q + E;    // total losses
  return(NumericVector::create(L,E,0.0,qcap,Q));
}

//Feedback model needs function Q_E (bottom flux at the root zone) #equation 6 VvdZ2009
//do we still need this, same as Fun_E[5]?
double Q_E_cpp(double s,double Z,List soilpar,List vegpar) {
  double Zr = vegpar["Zr"];
  double b = soilpar["b"];
  double Ks = soilpar["K_s"];
  double hb = soilpar["hb"];
  double H1 = G_cpp(b,hb,(Z-Zr));
  double Q = Ks*(pow(s,(2*b+3))- H1); // bottom flux
  return(Q);
}

// Deep root functions              #VvdZ 2009 eq.7
// in RootFunctions.cpp

// Feedback DR Model Vervoort and van der Zee (2009)
// [[Rcpp::export]]
NumericVector Fun_FB_cpp(double s,List soilpar,List vegpar,double Z,double Zmean) {
  double Zr = vegpar["Zr"];
  double c1 = vegpar["c1"]; 
  double fs = vegpar["fs"]; 
  
  // Calculate the Transpiration components
  double E_max = E_Teuling_cpp(1,vegpar);
  double Rc = Rc_B_cpp(Z,Zr,c1,Zmean,fs);
  //    double fr = int_test(Zr,0.0,c1);
  double T_Zr = int_test(c1,Z,Zr)*E_Teuling_cpp(s,vegpar); 
  double T_DR = std::min(Rc*E_max,E_max-T_Zr); //eq3 Verv.2009
  // Qbottom part
  double k1 = Q_E_cpp(s,Z,soilpar,vegpar);    
  //k is the bottom flux at root zone. Positive means leakage. Negative k means capillary rise
  double qcap = 0.0;
  if (k1 < 0) qcap = abs(k1);
  double L = k1+T_Zr; // Loss of soil moisture. k is leakage, T_Zr is uptake by plants.
  NumericVector out = NumericVector::create(L,T_Zr,T_DR,qcap,k1); // returns the total losses, the T loss from soil (T_Zr) and the T from groundwater (T_DR)
  return out;
}                        

//-----------------------------------------------
// New function including deep roots and feedback
// added 13/09/2012
// Combines old rho_new1 with FB_fun
// Including Joep's adjustment for Z.prev
// [[Rcpp::export]]
NumericVector FB_new_cpp(double s,List soilpar_in,List vegpar_in, double Zfb, 
                         double Zmean, double Z_prev) {
  double Zr = vegpar_in["Zr"];
  double ss = vegpar_in["s_star"];
  double sw = vegpar_in["s_w"]; 
  double c1 = vegpar_in["c1"]; 
  double fs = vegpar_in["fs"]; 
  double n = soilpar_in["n"];
  double b = soilpar_in["b"];
  double hb = soilpar_in["hb"];
  double beta = soilpar_in["beta"]; 
  double s_lim = pow(((Zfb - Zr)/hb),(-1/b));
  // apply G function
  double Gfun = G_cpp(b,hb,(Zfb-Zr));
  // calculate parameters
  NumericVector m_values = m_fun(vegpar_in, soilpar_in, Zfb, Gfun);
  const double E_max = E_Teuling_cpp(1,vegpar_in);
  const double eta = E_max/(n*Zr);
  double Q = 0.0;
  double E = 0.0;
  // deep root definitions
  double Rc = Rc_B_cpp(Zfb,Zr,c1,Zmean,fs);
  // define s_cr: s-critical
  double s_cr = m_values[2]/eta*(ss - sw) + sw;
  //Rcpp::Rcout << "m_values[3] " << m_values[2] << std::endl;
  //Rcpp::Rcout << "s_cr " << s_cr << std::endl;
  
  
  // Now calculate the soil moisture
  if (s > s_lim) {
    Q = m_values[0]*(std::exp(beta*(s - s_lim))-1);
    E = eta;
    //Rcpp::Rcout << s << std::endl;
    //Rcpp::Rcout << Q << std::endl;
  } else {
    if (m_values[1] < eta) {
      
      if (s < s_cr) {
        NumericVector m_values_prev = m_fun(vegpar_in, soilpar_in, Z_prev, G_cpp(b,hb,(Z_prev-Zr)));
        //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
        E = 0.0;  
        Q = -(m_values[2]*((s_cr-sw)/(ss-sw))-m_values_prev[2]*((s-sw)/(ss-sw)));
        //Rcpp::Rcout << s << std::endl;
        //Rcpp::Rcout << Q << std::endl;
      }
      if (s > s_cr && s <= ss)
      {
        E = eta*((s-s_cr)/(ss-s_cr));
        Q = -m_values[2]*((s-s_cr)/(ss - s_cr));
        //Rcpp::Rcout << s << std::endl;
        //Rcpp::Rcout << Q << std::endl;
      }
      if (s > ss && s <= s_lim) {
        Q = -m_values[1]*(1-std::exp(beta*(s-s_lim)));
        E = eta;
        //Rcpp::Rcout << s << std::endl;
        //Rcpp::Rcout << Q << std::endl;
      }
    } else {
      if (s > ss && s <= s_lim) {
        E = eta;
        Q = -m_values[1]*exp(beta*(s-s_lim));
        //Rcpp::Rcout << s << std::endl;
        //Rcpp::Rcout << Q << std::endl;
      } else {    // again, if GW table rises: E=0 and qcap =difference qcap and old transpiration
        NumericVector m_values_prev = m_fun(vegpar_in, soilpar_in, Z_prev, G_cpp(b,hb,(Z_prev-Zr)));
        //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
        E = 0.0;  
        Q = -(m_values[2]*((s_cr-sw)/(ss-sw))-m_values_prev[2]*((s-sw)/(ss-sw)));
        // Rcpp::Rcout << s << std::endl;
        //Rcpp::Rcout << Q<< std::endl;
      }
    }
  }
  double T_Zr = int_test(c1,Zr,0)*E*(n*Zr); 
  double T_DR = std::min(Rc*E_max,E_max-T_Zr); //eq3 Verv.2009
  double Q_out = Q*n*Zr;   
  double L = Q_out + T_Zr;
  double qcap = 0.0;
  if(Q_out < 0.0) qcap = abs(Q_out);
  NumericVector out = NumericVector::create(L,T_Zr,T_DR,qcap,Q_out); 
  return out;
}

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
                      std::string model,
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


// /*** R
// timesTwo(42)
// */