# include <Rcpp.h>
// this does not work,check out:
// http://stackoverflow.com/questions/13995266/using-3rd-party-header-files-with-rcpp
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
                   double deltat,
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
    double ss_t = last_t_soil_sat + Phi - loss[0]/(n*Zr)*deltat;
    
// build in safety if capillary flux fills soil with dynamic water table
    if (ss_t > 1) {  // very wet conditions
      soil_sat = 1;
      qcap = (loss[3] - (ss_t - 1)*n*Zr/deltat)*deltat;
    } else { // normal conditions
      soil_sat = ss_t;
      qcap = loss[3]*deltat;
    }
    static_stress = zeta(soil_sat, ss, sw, q);
    Tg = loss[2]*deltat;
// Rcpp::Rcout << loss[2];
      Ts = loss[1]*deltat;
      Tt = Tg + Ts;
      if (loss[4] > 0.0) {
        Leakage = loss[4]*deltat;
      } 
      //increase with difference between originally calculated and new (pot.reduced) capillary flux
      smloss = (loss[0]*deltat  - (qcap - loss[3]*deltat)); 
      GWrech = Leakage-loss[2]*deltat-qcap;
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
NumericMatrix WBEcoHyd(NumericVector R, NumericVector ETp,
              CharacterVector vtype, List soilpar,
              NumericVector s_init,
              double fullday, 
              //std::string model,
              NumericVector Zmean,
              NumericVector GWdepths,
              NumericVector GWdepths_prev,              
              double deltat,
              NumericVector NX, NumericVector NY,
              int t) {
  // both R and Etp are single vectors of the rainfall and the ETpotential
  // vtype is a vector of vegetation types along the grid
  // soilpar is the specification of the soiltype (currently 1 type across the grid)
  // s_init are the inital soil moisture levels for the day
  // deltat is the increment in the integration

  // Zmean is a vector of long term mean gw depths
  // GWdepths is the gw depths for each gridcell at t
  // GWdepths_prev are the previous gw depths for each gridcell (at t-1)
  
  // define variables
  int n_int = 1/deltat;
  int m = NX.size()*NY.size();
  
  double R_in = 0.0;
 
  // storage output vectors
  // grid vectors
  NumericVector s_g(m), qcap_g(m), Tgw_g(m), Tsoil_g(m), Ttotal_g(m), Leakage_g(m);
  NumericVector GWrech_g(m), smloss_g(m), static_stress_g(m), surfoff_g(m);
  //Matrices
  NumericMatrix Storage_Subday(12,10);
  // this only works as long as NX = 1
  NumericMatrix Storage_Grid(10,NY.size());
  
  // run through grid
  for (int j = 0;j < m;j++) {
    // Call the vegpar function
    List vegpar = Veg_cpp(vtype[j], soilpar);
    double Ep = ETp[j];
    

    // define soil and vegparameters
    double Zr = vegpar["Zr"];
    double n = soilpar["n"];
    
    // Now run integration loop
    for (int p = 0;p < n_int;p++) {
      if (p == 0) {
        double R_in = R[t]/(n*Zr);
      } 
      
      // run the water balance
        double s_old = s_init[j];
        if (p > 1 ) s_old = Storage_Subday((p-1),1);
        // Storage_Subday is not defined
        List Water_out = WB_fun_cpp(vegpar, R_in, s_old, 
                 soilpar, Zmean[j], deltat, GWdepths_prev, GWdepths);
//Write water balance list output to subdaily vector	
        int xsize =  Storage_Subday.ncol();
        for (int k = 0; k < xsize; k++) {
          Storage_Subday(p,k) = Water_out[k];
        }
        // 
// store s_init for gridcell till next day
        if (p == 12) s_init(j) = Storage_Subday(p,1);
    } // close p loop
// recalculate to daily output and write away
// produce different output, produce daily values for each grid cell
      for (int k = 0; k < Storage_Grid.ncol(); k++) {
        Storage_Grid(j,k) = sum(Storage_Subday(_,k));
        if (k > 0) Storage_Grid(j,k) = mean(Storage_Subday(_,k));
      }
    //Rcpp::Rcout << "j =" << j;
  } // close j loop
  return Storage_Grid;
}