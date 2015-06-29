# Script to massage river data for Yule river
# change raw to daily and monthly data
# Ecohydrology 2 D modelling Elise Honours project
# 2012 09 28

# grab the current date and store
today <- format(Sys.Date(),"%Y%m%d")

# set working dir
setwd("f:/willem/research/ecohydrology/2dmodelling")

# read in the data, but skip 10 lines for header and strange data
Raw <- read.csv("data for elise/20120920_StreamHeightYuleRiver.csv",skip = 10,
	header=F)
# big file, will read slow
# check top of file
head(Raw)
# add columnnames
colnames(Raw)<- c("DateTime","Height","Quality")

# Create a column with just the date and no times
Raw$Date <- strptime(Raw$DateTime,"%d/%m/%Y")
# add a column of months
Raw$Month <- as.numeric(substr(Raw$Date,6,7))

# Now use aggregate to summarise by day
YuleDay <- aggregate(Raw$Height,list(Date=as.Date(Raw$Date)),mean,na.rm=T)
# quick plot. seems that the null level was adjusted several times
plot(YuleDay,type="l",xlab="Date",ylab="Stage height (m)")

# we only need the last couple of years
Yuleday.short <- YuleDay[YuleDay$Date > "2005-01-01",]
# plot and add columnnames
plot(Yuleday.short,type="l",xlab="Date",ylab="Stage height (m)")
head(Yuleday.short)
colnames(Yuleday.short) <- c("Date","Height")

# Now substract the minimum to zero the series
Yuleday.short$Height <- Yuleday.short$Height - min(Yuleday.short$Height,na.rm=T)
plot(Yuleday.short,type="l",xlab="Date",ylab="Stage height (m)")

# write to a file
write.csv(Yuleday.short,paste(today,"YuleRivLev.csv",sep="_"),row.names=F)