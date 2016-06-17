# Do some timeseries analysis as suggested by Toon and Sjoerd 20/01/2014

# 2-D modelling Ecohydrology
require(lattice)

setwd("x:/vervoort/research/ecohydrology/2dmodelling")
today <- format(Sys.Date(),"%Y%m%d")


#### Read in the Rain data and climate data
# adjusted WV 20120830
Climate <- read.csv("YuleClim.csv")
Rain <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), Rain=Climate[,2]/10) # already in cm/d
ETp <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), ETp=rep(Climate[,5]/10))
GWdata <- read.csv("20131213_YuleRiverGWdata.csv")
GWdata.adj = data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"),
                        GWlevel=as.numeric(GWdata[match(as.Date(Climate[,1], "%d/%m/%Y"),
                                                   as.Date(GWdata[,1], "%d/%m/%Y")),2]))
# adjusted WV 20120928 
# daily streamflow
Stream <- read.csv("20120928_yulerivlev.csv")
# This is the righthand boundary condition

# There are missing data in the stream data
# need to match the stream dates to the climate dates
Stream.adj <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"),
                         Height=Stream$Height[match(as.Date(Climate[,1], "%d/%m/%Y"),
                                                    as.Date(Stream[,1], "%d/%m/%Y"))])

# correlations, autocorrelations and cross correlations

# rainfall
acf(Rain$Rain)
# Somewhat correlated, periodic, as you would expect from the monsoonal characteristic
# groundwater, is this equally spaced?
acf(GWdata.adj$GWlevel,na.action=na.pass)
# not very correlated, purely due to missing data
acf(Stream.adj$Height,na.action=na.pass)
# highly correlated as would be expected

# correlations
plot(Rain$Rain,GWdata.adj$GWlevel)
# no relationship
plot(Stream.adj$Height,GWdata.adj$GWlevel)
# some relationship as expected

# cross correlations
ccf(GWdata.adj$GWlevel,Rain$Rain,na.action=na.pass)
# rain leads gw at short lags and between 15 and 20 day lags
ccf(GWdata.adj$GWlevel,Stream.adj$Height,na.action=na.pass)
# much stronger cross correlation, again at short lags and at 10 to 20 days

