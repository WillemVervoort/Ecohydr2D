# generate all input data
require(lattice)
setwd("x:/vervoort/research/ecohydrology/2dmodelling")

#### Read in the Rain data and climate data
Climate <- read.csv("2dModellingBaseDataYuleClim.csv")
Rain <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), Rain=Climate[,2]/10) # already in cm/d
# Put in no rain
#Rain <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), Rain=rep(0,nrow(Climate))) 
ETp <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), ETp=rep(Climate[,5]/10))
GWdata <- read.csv("2dModellingBaseData/YuleRiverGWdata.csv")

#xyplot(Rain~dates,data=Rain,type="h")

# daily streamflow
Stream <- read.csv("2dModellingBaseData/yulerivlev.csv")
# This is the righthand boundary condition
head(Stream)
xyplot(Stream$Height~as.Date(Stream$Date,"%d/%m/%Y"),type="l")

# There are missing data in the stream data
# need to match the stream dates to the climate dates
Stream.adj <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"),
                         Height=Stream$Height[match(as.Date(Climate[,1], "%d/%m/%Y"),
                                                    as.Date(Stream[,1], "%d/%m/%Y"))])

# now need to fill the missing values as the model does not deal with NA
# use spline.fun with rainfall to fill
st.fun <- splinefun(Climate[,2],Stream.adj[,2])
Stream.adj[is.na(Stream.adj[,2]==T),2] <- st.fun(Climate[is.na(Stream.adj[,2]==T),2])

z <- nchar(as.character(Stream[,1]))
Stream.month <- aggregate(Stream$Height,
                          list(month=substr(Stream[,1],z-6,z-5),
                               year=substr(Stream[,1],z-3,z)),
                          sum)
head(Stream.month)
Stream.month$date <- as.Date(paste(Stream.month$year,Stream.month$month,"01",
                                   sep="-"))
xyplot(x~date,Stream.month,type="h",lwd=3,col="gray25",xlab=list("Date",font=2,cex=1.2),
       ylab=list("Monthly streamflow in mm",font=2,cex=1.2),scales=list(font=2,cex=1.2))

save.image("2dModellingBaseData/Allinput.rdata")
