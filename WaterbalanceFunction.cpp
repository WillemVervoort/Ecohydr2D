# include <Rcpp.h>
using namespace Rcpp;

// Rewriting the 2-D water balance implementation across grid cells
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

// E_Teuling
double E_Teuling(double s, List vegpar) {
	if(s <= vegpar["s_w"]) {
		double beta_T = 0.0;
	} else {
		if (s <= vegpar["s_star"] & s > vegpar["s_w"]) {
			beta_T = (s - vegpar["s_w"])/(vegpar["s_star"] - vegpar["s_w"]);
		} else {
			beta_T = 1.0;
		}
	}
	double E = vegpar["fr"]*beta_T*(1 - exp(-vegpar["c_T"]*vegpar["LAI"]))*vegpar["Ep"];
	return(E);
}	
	
// G function
double G(double b, double hb, double Z) {
	double b1 = 2.0 + 3.0/b;
	double a1 = 1.0 + (3.0/2.0)/(b1 - 1.0);
	double H1 = a1*pow(hb/Z,b1);
}

// rho_new_1
NumericVector rho_new_1(double s, double Z, List soilpar, List vegpar,
			Z.mean=NULL, Z.prev=NULL) {
				// define variables
				double Zr = vegpar["Zr"];
				double n = soilpar["n"]
				double b = soilpar["b"]
				double Ks = soilpar["K_s"]
				double hb = soilpar["psi_s_bar"]*(-10E4);
				double s.lim = pow((Z - Zr)/hb,(-1/b));
				// apply G function
				double G1 = G(b,hb,(Z-Zr));
				// calculate parameters
				double m = Ks/(n*Zr*(exp(soilpar["beta"]*(1-s.lim))-1));
				double m1 = Ks*G1/(n*Zr*(exp(soilpar["beta"]*(vegpar["s_star"]-s.lim))-1));
				double m1 = Ks*G1/(n*Zr);
				double E_max = (1-exp(-vegpar["c_T"]*vegpar["LAI"]))*vegapr["Ep"];
				double eta = E_max/(n*Zr);
				double r = 0;
				
				// Now calculate the soil moisture
				if (s > s.lim) {
					double Q = m*(exp(soilpar["beta"]*(s - s.lim))-1);
					double E = eta;
				}
				if (m1 < eta) {
					// define s_cr: s-critical
					double s_cr = m2/eta*(vegpar["s_star"] - vegpar["s_w"]) + vegpar["s_w"];
					
					
				}
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