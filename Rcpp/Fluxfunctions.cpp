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
    return std::exp(-c1*z/100);
  }
};

// equation 8 VvdZ 2009
// [[Rcpp::export]]
double int_test(double c1,double z1,double z2) {
  test_f f(c1);
  double err_est;
  int err_code;
  double res = integrate(f,z2,z1, err_est, err_code)/integrate(f,0,1e5, err_est, err_code);
  return res;
}

// new deep root functions following Orellana et al. 2012
double s_fun_cpp(double x,double fs) {
  return exp(-(fs*x)); 
}

// Zmean is the average long term groundwater depth
double  RWU_cpp(double z1,double Zmean, double fs) {
  return R::dnorm(z1/100,Zmean/100,s_fun_cpp(Zmean/100,fs),FALSE)*2; // test values
} 
// [[Rcpp::export]]
double Rc_B_cpp(double z1,double z2,double c1,double Zmean,double fs) {
  return RWU_cpp(z1,Zmean,fs)*int_test(c1,z1,z2);
}
//--------------------------------------------------
#include <Rcpp.h>
using namespace Rcpp;

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
NumericVector FB_new_cpp(float s,List soilpar_in,List vegpar_in, double Zfb, 
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

/*** R
sourceCpp("soilfun.cpp")
sourceCpp("Vegfun.cpp")
sourceCpp("RootFunctions.cpp")
soilpar = Soil_cpp("M Clay")
vegpar = Veg_cpp("TreesDR", soilpar)
FB_new_cpp(0.4,soilpar,vegpar,200,190,190)
*/
