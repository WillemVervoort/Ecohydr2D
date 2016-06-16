setwd("x:/vervoort/research/rcode/ecohydrology/ecohydrology2dmodellng")

require(Rcpp)
sourceCpp("rcpp/test.cpp")

#sourceCpp("rcpp/soilfun.cpp")

#sourceCpp("rcpp/Vegfun.cpp")

sourceCpp("rcpp/Rootfunctions.cpp")

sourceCpp("rcpp/Fluxfunctions.cpp")
