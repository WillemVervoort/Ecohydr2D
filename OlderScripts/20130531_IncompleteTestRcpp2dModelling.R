# rcpp script for Joep's function

# ---define and compile the c++ code--------------
library(inline)
src <- '

// inputs are: N,stype,vtype,Rain,ETp,stream,gw_in,Zmean,pumprate=0, NX, NY
  int N = Rcpp::as<int>(N_In);
  Rcpp::CharacterVector stype(stype_In);
  Rcpp::CharacterVector vtype(vtype_In);
  Rcpp::NumericMatrix stream(stream_In);
  Rcpp::NumericVector gw_in(gw_In);
  double Zmean = Rcpp::as<double>(Zmean_In);
  double pumprate = Rcpp::as<double>(pumprate_In);
  int NX = Rcpp::as<int>(NX_In);
  int NY = Rcpp::as<int>(NX_In);
  double deltat = 1./12.;
  int Rowsize = as.int(N * 1./deltat);
//Internal vectors
  Rcpp::NumericMatrix ss(Rowsize, NX*NY);
  Rcpp::NumericMatrix Phi(Rowsize, NX*NY);
  Rcpp::NumericMatrix Tg(Rowsize, NX*NY); // Transpiration from deep roots
  Rcpp::NumericMatrix qcap(Rowsize, NX*NY);
  Rcpp::NumericMatrix Ts(Rowsize, NX*NY);
  Rcpp::NumericMatrix Tt(Rowsize, NX*NY);
  Rcpp::NumericMatrix Leakage(Rowsize, NX*NY);
  Rcpp::NumericMatrix smloss(Rowsize, NX*NY);
  Rcpp::NumericMatrix surfoff(Rowsize, NX*NY);
// some extras to store internal data
  double t = NumericVector(Rowsize) // time vector
  int finalN = N //this will be overwritten when formula stops before N is reached  
  Rcpp::NumericVector loss_an(5);


 
// old rubbish
//  if (DR == TRUE) model <- "FB_new" else model <- rho_new_1 # "Fun_E" #
//  sinit<-vegpar$s_star
//  Zr <- vegpar$Zr
//  n <- soilpar$n
//  spec_y <- soilpar$spec_y
//  R <- Rain[,2] #vector with rain data
//  DELTvector<-R   # vector with DELT values to see when it is calculated again. Same length as Rain


//double n = soilpar(0);
//double  Zr = vegpar(3);
//int s_size = ss.size(); // calculates the size of ss??
// initialise soil saturation
ss(0) = 0.5;

for(int i = 0; i < (s_size-1); i++){
     // make sure all rainfall occurs in first integration increment
	double first = floor(float(i)/float(deltat));
     	if(first == (float(i)/float(deltat))){
		if((R(first)/(n * Zr)) > (1.0 - ss(i))) {
			// Calculate infiltration
			Phi(i) = (1.0 -ss(i));
			// calculate runoff
			 answer(i,7) = R(first)-(n * float(Zr)) * (1. - ss(i));
		} else {
			 Phi(i) = R(first)/(n * float(Zr));
      	     answer(i,7) = 0.0;
        	}
	} else {
		Phi(i) = 0.0;
	      answer(i,7) = 0.0;
	}
	// call loss model
	loss_an = loss(ss(i), soilpar,vegpar, Z, mpar);
	// calculate water balance
	ss(i+1) = ss(i) + Phi(i) - loss_an(0) / float(deltat) / (n * float(Zr));
	// Store results
	answer(i,0) = ss(i+1); // soil water
	answer(i,1) = loss_an(1)/float(deltat); // soil root zone T
	answer(i,2) = loss_an(2)/float(deltat); // deep root T
    	answer(i,3) = loss_an(3)/float(deltat); // capillary flux
    	answer(i,4) = loss_an(0)/float(deltat)-loss_an(1)/float(deltat);// drainage
    	answer(i,5) = Phi(i)/float(deltat); // unadjusted infiltration
    	answer(i,6) = loss_an(4)/float(deltat); //unadjusted bottom q 

}
		
return wrap(answer);
'
# Included subroutines loss function and et function
inc <- '
	double et(double s, Rcpp::NumericVector vegpar){
// Teuling and Troch (2005) function
	double E, beta_T;
     // macroscopic loss function
	if( s <= vegpar(5) ) {
		beta_T = 0.0;
		} else {
			if(s <=vegpar(6) && s > vegpar(5)){
				beta_T = (s - vegpar(5))/(vegpar(6)-vegpar(5));
			 } else {
         		beta_T = 1.0;
			}
		}			
     // calculate actual ET
	E = vegpar(9)*beta_T*(1.0-exp(-vegpar(8) * vegpar(7))) * vegpar(10);
	return E;
}

// loss function
Rcpp::NumericVector loss(double s, Rcpp::NumericVector soilpar,Rcpp::NumericVector vegpar, double Z, Rcpp::NumericVector mpar){
	double hb,b1,a1,H1,Qb,E_max, E, T_Zr, T_DR;
	Rcpp::NumericVector loss_a(5);
      // calculate parameters campbell loss function
	 hb = -pow(10.0,4.0) * soilpar(3);
       b1 = 2.0 + 3.0/soilpar(2);
       a1 = 1.0 + (3.0/2.0)/(b1 - 1.0);
	// Eagleson function
       H1 = a1 * pow((hb/ float(Z - vegpar(3) ) ), b1);
	// Calculate max ET and actual ET
		E_max =(1.0-exp(-vegpar(8) * vegpar(7))) * vegpar(10);
		// call et function
		E =  et(s,vegpar);
		// transpiration from rootzone
		T_Zr = mpar(16)*E;
		loss_a(1) = T_Zr;
		// transpiration deep roots
		if(mpar(15) * E_max <= E_max - T_Zr){
			T_DR =	mpar(15) * E_max;
		} else {
			T_DR =	E_max - T_Zr;
		}
		loss_a(2) = T_DR;
		// bottom flux
		Qb = soilpar(1) * (pow(s,(2.*soilpar(2) + 3.)) - H1);
		loss_a(4) = Qb;
		// calculate capillary component
		if (Qb <= 0){
	         loss_a(3) = fabs(Qb);
      	 } else {
          	   loss_a(3) = 0.0;
		}
       	loss_a(0) = T_Zr + Qb;// total loss from rootzone
       
       return loss_a;
       }
       
'
cfun <- cxxfunction(signature(R_In="numeric",Z_In="int",
            time_In="any",deltat_In="any",
            soilpar_In="any",vegpar_In="numeric",mpar_In="numeric"), 
            src, plugin="Rcpp", includes=inc,verbose =F)
# ------------------------------------------------ end compiling function
# all above can be moved to another script and called once via source()
# only works if R is in dir with no "space" in the dir name
# also only works if Rtools is installed

WBcxx <- function(R,Z,time,deltat,soilpar,vegpar,mpar) {
	run <- cfun(R_In=rf, Z_In=325, time_In=time,  
        deltat=deltat, soilpar_In=soilparVec, 
        vegpar_In=vegparVec,mpar_In=mparVec)
  return(data.frame(s=run[,1],TZr=run[,2],Tg=run[,3],
	qcap=run[,4],Drain=run[,5],Inf_u=run[,6],Bflux=run[,7],Q=run[,8])) 

}	


