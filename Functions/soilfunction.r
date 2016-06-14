#
# Soil data as a function
# Ecohydrology WIMEK project
# Willem Vervoort September 2007
# To simplify soil input
# ########################################
Soil <- function(stype) {
   psi_sh<--10
   soil <- stype
  if (soil == "L Med Clay") {
    # Medium Light Clay
    n<-0.418 # porosity
    # more soil variables for evaporation & losses
    K_s<-3.51 # cm/day
    b<-13.48 # neurotheta LMC
    nvg <- 1.089
    avg <- 0.0591
    s_fc<-0.364/n # Field capacity
    psi_s_bar<--1.5E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    # This is the bubbling pressure
    spec_y <- 0.054 #Johnson 1967 says 0.05, specific yield. 
  }
  if (soil == "L Med Clay Stony") {
    # Medium Light Clay
    n<-0.318 # porosity
    # more soil variables for evaporation & losses
    K_s<-3.51 # cm/day
    b<-13.48 # neurotheta LMC
    nvg <- 1.089
    avg <- 0.0591
    s_fc<-0.264/n # Field capacity
    psi_s_bar<--1.5E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    # This is the bubbling pressure
    spec_y <- 0.054 #Johnson 1967 says 0.05, specific yield. 
  }
                      
  if (soil == "S Clay Loam") {
    # Sandy Clay Loam
    n<-0.367 # porosity
    # more soil variables for evaporation & losses
    K_s<-52.08 # cm/day
    b<-6.4069 # neurotheta sandy clay loam
    avg <- 0.0521
    nvg <- 1.237
    s_fc<-0.2677/n # Field capacity

    psi_s_bar<--1.2E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    spec_y <- 0.07  #difference Fc and por Johnson 1967 says 0.07 
  }

  if (soil == "Loamy Sand") {
    # Loamy sand
    n<-0.37 # porosity
    # more soil variables for evaporation & losses
    K_s<-175.3 # cm/day
    b<-4.5206
    avg <- 0.0641
    nvg <- 1.344
    s_fc<-0.2098/n # Field capacity

    psi_s_bar<--0.66E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    # This is the bubbling pressure
    
    spec_y <- 0.17  # changed to 0.1 to increase rise in gw, not difference por and fc
  }

  if (soil == "H Clay") {
    # Medium Heavy Clay
    n<-0.4473 # porosity
    # more soil variables for evaporation & losses
    K_s<-2.82 # cm/day
    b<-16.1501 # neurotheta Medium heavy clay
    avg <- 0.0613
    nvg <- 1.086
    s_fc<-0.3936/n # Field capacity

    psi_s_bar<--1.4E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    spec_y <- 0.05  # difference por and fc
  }

  if (soil == "M Clay") {
    # Medium Clay
    n<-0.4391 # porosity
    # more soil variables for evaporation & losses
    K_s<-6.04 # cm/day
    b<- 13.5127
    avg <- 0.0507
    nvg <- 1.088
    s_fc<-0.3818/n # Field capacity

    psi_s_bar<--1.75E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    # This is the bubbling pressure
    spec_y <- 0.05 # difference por and fc
  }
 if (soil == "C Sand") {
    # Coarse Sand
    n<-0.368 # porosity
    # more soil variables for evaporation & losses
    K_s<-182.68 # cm/day
    b<- 4.1152
    avg <- 0.0712
    nvg <- 1.392
    s_fc<-0.1895/n # Field capacity

    psi_s_bar<--0.61E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    # This is the bubbling pressure
    spec_y <- 0.27  # difference por and fc
  }


# Other derived parameters
    s_h<-(psi_sh/psi_s_bar)^(-1/b)# soil moisture at hygroscopic point
    beta<-2*b+4 #page 714 Laio et al., 2001a

    # Define parameters for Eagleson function
    # Beta parameters
    beta1 <- 2+3/b
    # alpha parameter
    a1 <- 1+(3/2)/(beta1-1)
    soilpar <- list(n=n,K_s=K_s,b=b,psi_s_bar=psi_s_bar,s_fc=s_fc,s_h=s_h,beta=beta,beta1=beta1,a1=a1,nvg=nvg,avg=avg,spec_y=spec_y)
    return(soilpar)
}

