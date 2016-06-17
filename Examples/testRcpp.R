setwd("x:/vervoort/research/rcode/ecohydrology/ecohydrology2dmodellng")

require(Rcpp)
sourceCpp("rcpp/test.cpp")

#sourceCpp("rcpp/soilfun.cpp")

#sourceCpp("rcpp/Vegfun.cpp")

sourceCpp("rcpp/Rootfunctions.cpp")

sourceCpp("rcpp/Fluxfunctions.cpp")


# check against the R implementation
source("functions/soilfunction.r")
source("functions/vegfunction.r")
source("functions/FluxFunctions2DModel.R")

FB_new_cpp(0.4,soilpar,vegpar,200,190,190)
# R implementation
soilpar_r <- Soil("M Clay")
vegpar_r <- Veg("TreesDR", soilpar_r)
FB_new(0.4,soilpar_r,vegpar_r,200,190,190)
