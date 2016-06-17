# generate ET series for Garbiel for 5 years

# based on data series from Sebastian

setwd("x:/vervoort/research/ecohydrology/2dmodelling")

ETdata.Seb <- read.csv("20130624_TranspirationTrees.csv")
head(ETdata.Seb)

require(zoo)

ET.seb <- zoo(ETdata.Seb[,2],order.by=as.Date(ETdata.Seb$Date,"%d/%m/%Y"),frequency=1)

plot(ET.seb)
times <- julian(as.Date(ETdata.Seb$Date,"%d/%m/%Y"),origin=as.Date("2008-01-01"))

set.seed(4356)

# test function
cos.ser <- 2+0.3*cos(2*pi*times/365) + rnorm(times,0,0.15)

# real function
T.sim <- function(times) {
  return(2+0.3*cos(2*pi*times/365) + rnorm(times,0,0.15))
}

# generate 5 years of virtual transpiration data
# some dates
dates <- seq.Date(as.Date("2008-01-01"),to=as.Date("2013-12-31"),by=1)
# number of days
days <- 1:length(dates) 

# generate transpiration
Transp.s <- T.sim(days)

plot(dates,Transp.s,type="l",ylab="simulated transpiration (mm/day)",
     xlab="Dates")

# write away in a file
write.csv(data.frame(Date=dates,Transp=Transp.s),
          "20131105_SimulatedTranspiration5years.csv",row.names=F)
