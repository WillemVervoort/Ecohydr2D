## zzz.R
## Following example from hydromad: 
# Hydrological Modelling and Analysis of Data
## Willem Vervoort
##
##

## local environment to store user options
.ecohydro2dEnv <- new.env()
.ecohydro2dEnv$options <- list()

.onLoad <- function(libname, pkgname)
{
  ecohydro2d.options(.defaultecohydro2dOptions())
}