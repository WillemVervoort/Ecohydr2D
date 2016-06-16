// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;


// from rootfunctions.cpp to be removed when this is a package
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


double  RWU(double z1,double Zmean, double fs) {
  return R::dnorm(z1/100,Zmean/100,s_fun(Zmean/100,fs),FALSE)*2; // test values
} 
// [[Rcpp::export]]
double Rc_B_cpp(double z1,double z2,double c1,double Zmean,double fs) {
  return RWU(z1,Zmean,fs)*int_test(c1,z1,z2);
}
//--------------------------------------------------
#include <Rcpp.h>
using namespace Rcpp;

// E_Teuling
double E_Teuling(double s, List vegpar) {
  double sw = vegpar["s_w"]; 
  double ss = vegpar["s_star"];
  double beta_T = 0.0;
  if(s > sw) {
    if ((s < ss) & (s > sw)) {
      beta_T = (s - sw)/(ss - sw);
    } else {
      beta_T = 1.0;
    }
  }
  double fr = vegpar["fr"]; 
  double c_T = vegpar["c_T"]; 
  double LAI = vegpar["LAI"]; 
  double Ep = vegpar["Ep"]; 
  double E = fr*beta_T*(1 - exp(-c_T*LAI))*Ep;
  return(E);
}	

// G function
// [[Rcpp::export]]
double G(double b, double hb, double Z) {
  double b1 = 2.0 + 3.0/b;
  double a1 = 1.0 + (3.0/2.0)/(b1 - 1.0);
  double H1 = a1*pow(hb/Z,b1);
  return(H1);
}

// function to generate m values
// [[Rcpp::export]]
NumericVector m_fun(List vegpar, List soilpar, double Z_in, double G1) {
  double Zr = vegpar["Zr"]; 
  double b = soilpar["b"]; 
  double hb = soilpar["hb"]; 
  double Ks = soilpar["K_s"]; 
  double n = soilpar["n"]; 
  double beta = soilpar["beta"]; 
  double ss = vegpar["s_star"];
  
  double s_lim = pow((Z_in - Zr)/hb,(-1/b));
  double m = Ks/(n*Zr*(exp(beta*(1-s_lim))-1));
  double m1 = Ks*G1/(n*Zr*(exp(beta*(ss-s_lim))-1));
  double m2 = Ks*G1/(n*Zr);
  return(NumericVector::create(m,m1,m2));
}

// rho_new_1
// [[Rcpp::export]]
NumericVector rho_new_1(double s, double Z, List soilpar, List vegpar,
                        double Z_mean, double Z_prev) {
  // define variables
  double Zr = vegpar["Zr"];
  double ss = vegpar["s_star"];
  double sw = vegpar["s_w"]; 
  double n = soilpar["n"];
  double b = soilpar["b"];
  double hb = soilpar["hb"];
  double beta = soilpar["beta"]; 
  double s_lim = pow((Z - Zr)/hb,(-1/b));
  // apply G function
  double G1 = G(b,hb,(Z-Zr));
  // calculate parameters
  NumericVector m_values = m_fun(vegpar, soilpar, Z, G1);
  double E_max = E_Teuling(s=1,vegpar);
  double eta = E_max/(n*Zr);
  double Q = 0.0;
  double E = 0.0;
  
  // Now calculate the soil moisture
  if (s > s_lim) {
    Q = m_values[1]*(exp(beta*(s - s_lim))-1);
    E = eta;
  }
  if (m_values[2] < eta) {
    // define s_cr: s-critical
    double s_cr = m_values[3]/eta*(ss - sw) + sw;
    
    if (s < s_cr) {
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G1 = G(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      E = 0.0;  
      Q = -(m_values[3]*((s_cr-sw)/(ss-sw))-m_values_prev[3]*((s-sw)/(ss-sw)));
    }
    if (s > s_cr && s <= ss)
    {
      E = eta*((s-s_cr)/(ss-s_cr));
      Q = -m_values[3]*((s-s_cr)/(ss - s_cr));
    } else {
      if (s > ss && s <= s_lim) {
        Q = -m_values[2]*(1-exp(beta*(s-s_lim)));
        E = eta;
      }
    }
  } else {
    if (s > ss && s <= s_lim) {
      E = eta;
      Q = -m_values[2]*exp(beta*(s-s_lim));
    } else {    // again, if GW table rises: E=0 and qcap =difference qcap and old transpiration
      double s_cr = m_values[3]/eta*(ss - sw) + sw;
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G1 = G(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      E = 0.0;  
      Q = -(m_values[3]*((s_cr-sw)/(ss-sw))-m_values_prev[3]*((s-sw)/(ss-sw)));
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
NumericVector Fun_E(double s,double Z,List soilpar,List vegpar) {
  double Zr = vegpar["Zr"];
  double b = soilpar["b"];
  double Ks = soilpar["K_s"];
  double hb = soilpar["hb"];
  double H1 = G(b,hb,(Z-Zr));
  double Q = Ks*(pow(s,(2*b+3))- H1); // bottom flux
  double  E = E_Teuling(s,vegpar);
  double qcap = 0.0;
  if (Q < 0) qcap = abs(Q); // capillary flow component
  double L = Q + E;    // total losses
  return(NumericVector::create(L,E,0.0,qcap,Q));
}

// Feedback model needs function Q_E (bottom flux at the root zone) #equation 6 VvdZ2009
// do we still need this, same as Fun_E[5]
// double Q_E(double s,double Z,List soilpar,List vegpar) {
//   double Zr = vegpar["Zr"];
//   double b = soilpar["b"];
//   double Ks = soilpar["K_s"];
//   double psi_s_bar = soilpar["psi_s_bar"];
//   double hb = psi_s_bar*(-10E4);
//   double H1 = G(b,hb,(Z-Zr));
//   double Q = Ks*(pow(s,(2*b+3))- H1); // bottom flux
//   return(Q);
// }

// Deep root functions              #VvdZ 2009 eq.7
// in RootFunctions.cpp

// Feedback DR Model Vervoort and van der Zee (2009)
// [[Rcpp::export]]
NumericVector Fun_FB(double s,List soilpar,List vegpar,double Z,double Zmean) {
  double Zr = vegpar["Zr"];
  double c1 = vegpar["c1"]; 
  double fs = vegpar["fs"]; 

// Calculate the Transpiration components
    double E_max = E_Teuling(1,vegpar);
    double Rc = Rc_B_cpp(Z,c1,Zr,Zmean,fs);
    //    double fr = int_test(Zr,0.0,c1);
    double T_Zr = int_test(c1,Z,Zr)*E_Teuling(s,vegpar); 
    double T_DR = std::min(Rc*E_max,E_max-T_Zr); //eq3 Verv.2009
    // Qbottom part
    double k1 = Fun_E(s,Z,soilpar,vegpar)[5];    
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
NumericVector FB_new(double s,List soilpar,List vegpar, double Z, 
              double Zmean, double Z_prev) {
  double Zr = vegpar["Zr"];
  double ss = vegpar["s_star"];
  double sw = vegpar["s_w"]; 
  double c1 = vegpar["c1"]; 
  double fs = vegpar["fs"]; 
  double n = soilpar["n"];
  double b = soilpar["b"];
  double hb = soilpar["hb"];
  double beta = soilpar["beta"]; 
  double s_lim = pow((Z - Zr)/hb,(-1/b));
  // apply G function
  double G1 = G(b,hb,(Z-Zr));
  // calculate parameters
  NumericVector m_values = m_fun(vegpar, soilpar, Z, G1);
  double E_max = E_Teuling(s=1,vegpar);
  double eta = E_max/(n*Zr);
  double Q = 0.0;
  double E = 0.0;
  // deep root definitions
  double Rc = Rc_B_cpp(Z,c1,Zr,Zmean,fs);
  // Now calculate the soil moisture
  if (s > s_lim) {
    Q = m_values[1]*(exp(beta*(s - s_lim))-1);
    E = eta;
  }
  if (m_values[2] < eta) {
    // define s_cr: s-critical
    double s_cr = m_values[3]/eta*(ss - sw) + sw;
    
    if (s < s_cr) {
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G1 = G(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      E = 0.0;  
      Q = -(m_values[3]*((s_cr-sw)/(ss-sw))-m_values_prev[3]*((s-sw)/(ss-sw)));
    }
    if (s > s_cr && s <= ss)
    {
      E = eta*((s-s_cr)/(ss-s_cr));
      Q = -m_values[3]*((s-s_cr)/(ss - s_cr));
    } else {
      if (s > ss && s <= s_lim) {
        Q = -m_values[2]*(1-exp(beta*(s-s_lim)));
        E = eta;
      }
    }
  } else {
    if (s > ss && s <= s_lim) {
      E = eta;
      Q = -m_values[2]*exp(beta*(s-s_lim));
    } else {    // again, if GW table rises: E=0 and qcap =difference qcap and old transpiration
      double s_cr = m_values[3]/eta*(ss - sw) + sw;
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G1 = G(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      E = 0.0;  
      Q = -(m_values[3]*((s_cr-sw)/(ss-sw))-m_values_prev[3]*((s-sw)/(ss-sw)));
    }
  }
  double T_Zr = int_test(c1,Z,Zr)*E*(n*Zr); 
  double T_DR = std::min(Rc*E_max,E_max-T_Zr); //eq3 Verv.2009
  double Q_out = Q*(n*Zr);   
  double E_out = E*(n*Zr);   
  double L = Q_out + E_out;
  double qcap = 0.0;
  if(Q_out < 0) qcap = abs(Q_out);
  NumericVector out = NumericVector::create(L,T_Zr,T_DR,qcap,Q_out); 
  return out;
}

// /*** R
// sourceCpp("soilfun.cpp")
// sourceCpp("Vegfun.cpp")
// sourceCpp("RootFunctions.cpp")
// soilpar = Soil("S Clay Loam")
// vegpar = Veg("TreesDR", soilpar)
// FB_new(0.4,soilpar,vegpar,200,190,190)
// */
