# include <Rcpp.h>
using namespace Rcpp;

// Rewriting the 2-D water balance implementation across grid cells

// E_Teuling
// [[Rcpp::export]]
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
  double hb = soilpar["hb"]; 
  double b = soilpar["b"]; 
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
NumericVector rho_new_1(double s, double Z, List soilpar, List vegpar,
                        double Z_mean=NULL, double Z_prev=NULL) {
  // define variables
  double Zr = vegpar["Zr"];
  double ss = vegpar["s_star"];
  double sw = vegpar["s_w"]; 
  double n = soilpar["n"];
  double b = soilpar["b"];
  double Ks = soilpar["K_s"];
  double psi_s_bar = soilpar["psi_s_bar"];
  double hb = psi_s_bar*(-10E4);
  double beta = soilpar["beta"]; 
  double s_lim = pow((Z - Zr)/hb,(-1/b));
  // apply G function
  double G1 = G(b,hb,(Z-Zr));
  // calculate parameters
  NumericVector m_values = m_fun(vegpar, soilpar, Z, G1);
  double E_max = E_Teuling(s=1,vegpar);
  double eta = E_max/(n*Zr);
  double r = 0.0;
  double Q = 0.0;
  double E = 0.0;
  
  // Now calculate the soil moisture
  if (s > s_lim) {
    double Q = m_values[1]*(exp(beta*(s - s_lim))-1);
    double E = eta;
  }
  if (m_values[2] < eta) {
    // define s_cr: s-critical
    double s_cr = m_values[3]/eta*(ss - sw) + sw;
    
    if (s < s_cr) {
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G1 = G(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      double E = 0.0;  
      double Q = -(m_values[3]*((s_cr-sw)/(ss-sw))-m_values_prev[3]*((s-sw)/(ss-sw)));
    }
    if (s > s_cr && s <= ss)
    {
      double E = eta*((s-s_cr)/(ss-s_cr));
      double Q = -m_values[3]*((s-s_cr)/(ss - s_cr));
    } else {
      if (s > ss && s <= s_lim) {
        double Q = -m_values[2]*(1-exp(beta*(s-s_lim)));
        double E = eta;
      }
    }
  } else {
    if (s > ss && s <= s_lim) {
      double E = eta;
      double Q = -m_values[2]*exp(beta*(s-s_lim));
    } else {    // again, if GW table rises: E=0 and qcap =difference qcap and old transpiration
      double s_cr = m_values[3]/eta*(ss - sw) + sw;
      NumericVector m_values_prev = m_fun(vegpar, soilpar, Z_prev, G1 = G(b,hb,(Z_prev-Zr)));
      //double s_cr_prev = m_values_prev[3]/eta*(ss - sw) + sw;
      double E = 0.0;  
      double Q = -(m_values[3]*((s_cr-sw)/(ss-sw))-m_values_prev[3]*((s-sw)/(ss-sw)));
    }
  }
    double Q_out = Q*(n*Zr);   
    double E_out = E*(n*Zr);   
    double L_out = Q_out + E_out;
    double qcap = 0.0;
    if(Q_out < 0) double qcap = abs(Q);
    return(NumericVector::create(L_out,E_out,0.0,qcap,Q_out));
}


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