# -------------------------------
# Some real scenario runs
# 
# This is the base case without rainfall as suggested by Sjoerd
# e-mail 20140109
require(lattice)

setwd("x:/vervoort/research/ecohydrology/2dmodelling")
today <- format(Sys.Date(),"%Y%m%d")


# In advance, define all parameters in gwt_input_parameters.r
#edit(file="X:/vervoort/research/rcode/ecohydrology/2dmodelling/gwt_input_parameters.r")


#### Read in the Rain data and climate data
# adjusted WV 20120830
Climate <- read.csv("YuleClim.csv")
Rain <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), Rain=Climate[,2]/10) # already in cm/d
# Put in no rain
#Rain <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), Rain=rep(0,nrow(Climate))) 
ETp <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), ETp=rep(Climate[,5]/10))
GWdata <- read.csv("20131213_YuleRiverGWdata.csv")

#xyplot(Rain~dates,data=Rain,type="h")

# adjusted WV 20120928 
# daily streamflow
Stream <- read.csv("20120928_yulerivlev.csv")
# This is the righthand boundary condition
head(Stream)
#xyplot(Stream$Height~as.Date(Stream$Date,"%d/%m/%Y"),type="l")

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

# Now run a few different simulations
# define initial gwheads and Zmean
#load("finalgwlevel.rdata")
#gw_in <- c(finalgw[1:(NX-1)],-12)
gw_in <-  GWdata[1,2] + 0.5 - 0.048*(1:42)#seq(-8,-10,length=42) #finalgw 
Zmean <- rep(600,42) #gw_in*-100 # could also try Zmean <- 500
fs_veg <- 0.25 # this is the spreading of the root water uptake around z_mean

# TO MODEL DEEP ROOTS WITH TREES YOU NEED TO CHANGE vtype to
# vtype should now be a vector
# vtype = "TreesDR" otherwise "TreesNoDR"
set.seed(1000)
veg <- rep("Grass",32)
veg[sample(1:length(veg),size=ceiling(0.3*length(veg)))] <- 
			"Bare"
veggies <- c(rep("TreesDR",10), veg)
soils <- "L Med Clay"
# define the soil for the aquifer
sp <- Soil_cpp("Loamy Sand") #Soil("L Med Clay")
#
# now includes separate K for aquifer (sand) based on Pfautsch et al.
result <- big.fun(N=nrow(Rain),stype=soils, aqK_in = sp$K_s/100, 
                  aq_specy_in  = sp$spec_y,vtype=veggies,Rain=Rain, 
                  ETp=ETp,stream=Stream.adj, gwheads = gw_in, Zmean = Zmean, 
                  today.m = today, fs = fs_veg)



# run stackfun on the output file
Out <- stackfun(result,42)
head(Out)

# put in dates rather than sequential numbers:
Out$Dates <- rep(Rain$dates,42)

# Now use Lattice for some pretty plotting
#Define a strip function
my.strip <- function(which.given, ..., factor.levels, bg) {
            levs <- paste(distancetoriver,"m")
            my.bg <- "gray90"
          strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}

# Set some background levels
trellis.par.set("regions", list(col=gray(25:100/100)))
trellis.par.set("background", list(col="white"))
# Make a histogram of the s values
#histogram(~s|as.factor(loc),data=Out,strip=my.strip, #subset=loc<100,
#	xlab="Soil saturation",ylab="Frequency")
#histogram(~qcap|as.factor(loc),data=Out,strip=my.strip, #subset=loc<100,
#	xlab="Capillary flow",ylab="Frequency")
#histogram(~gwrech|as.factor(loc),data=Out,strip=my.strip, #subset=loc<100,
#	xlab="groundwater recharge",ylab="Frequency")



# Make a histogram of Total Transpiration to compare to observed
#histogram(~(Tg*10)|as.factor(loc),data=Out,strip=my.strip,subset=loc<95,
#	xlab="GW transpiration (mm/day)",ylab="Frequency")

# read in the observed Transpiration
Tobs <- read.csv("20130624_TranspirationTrees.csv")



my.strip3 <- function(which.given, ..., factor.levels, bg) {
            levs <- paste("from river",c(0, 45, 85, 162.5),"m")
            my.bg <- "gray90"
          strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}


sub1_Out <- Out[(Out$loc %in% c(0, 45, 85, 162.5)),]

xyplot(Tg*10+Ts*10~as.Date(Dates)|as.factor(loc),data=sub1_Out, 
	type="l",strip=my.strip3,xlab =list( "Days",font=2,cex=1.2), 
	ylab = list("Transpiration (m)",font=2,cex=1.2),
	par.strip.text=list(font=2,cex=1.2), lty=c(1,2), col=c(1,"Gray50"),
	scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),as.table=T,
	key=list(text=list(lab=c("Groundwater","Soil"),cex=1),
        lines=list(lwd=2,lty=c(1,2),col=c(1,"Gray50")),x = .53, y = .35, corner = c(0, 0)))


GWsim <- data.frame(Dates=sub1_Out$Dates, level=sub1_Out$gwlevel,
                    ind=rep("sim",nrow(sub1_Out)),loc=sub1_Out$loc)
GWobs <- data.frame(Dates=rep(as.Date(GWdata$Date,"%d/%m/%Y"),4), 
                    level=rep(GWdata$Control_m_f_surf,4),
                    ind=rep("obs",nrow(GWdata)*4),
                    loc=rep(c(0, 45, 85, 162.5),each=nrow(GWdata)))
Strobs <- data.frame(Dates=rep(as.Date(Stream$Date,"%d/%m/%Y"),4), 
                    level=rep(Stream$Height,4),
                    ind=rep("stream",nrow(Stream)*4),
                    loc=rep(c(0, 45, 85, 162.5),each=nrow(Stream)))
# Rainobs <- data.frame(Dates=rep(as.Date(Rain$dates,"%d/%m/%Y"),4), 
#                      level=rep(Rain$Rain/2-5,4),
#                      ind=rep("rain",nrow(Rain)*4),
#                      loc=rep(c(0, 45, 85, 162.5),each=nrow(Rain)))

#GWplot <- rbind(GWsim,GWobs,Strobs,Rainobs)
GWplot <- rbind(GWsim,GWobs,Strobs)


xyplot(level~as.Date(Dates)|as.factor(loc),data=GWplot, groups=ind,
	type="l",strip=my.strip3,xlab =list( "Days",font=2,cex=1.2), 
	ylab = list("height",font=2,cex=1.2), col=c(1:3,"gray50"),
	par.strip.text=list(font=2,cex=1.2),#ylim=c(-11,-3),
	scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),as.table=T)

qqplot(Tobs$Control_Transpiration,Out[Out$loc==45,"Ttotal"]*10,
	xlab="Observed transpiration (sapflow)",
	ylab="Predicted total transpiration",
cex.lab=1.2,cex.axis=1.2,font.axis=2,font.lab=2)
lines(Tobs$Control_Transpiration,Tobs$Control_Transpiration,lty=2,lwd=2)

OutTt <- Out[Out$loc==0,"Ttotal"]
OutTt.red <- OutTt[OutTt > 0.1]

qqplot(Tobs$Control_Transpiration,OutTt.red*10,
	xlab="Observed transpiration (sapflow)",
	ylab="Predicted total transpiration",
cex.lab=1.2,cex.axis=1.2,font.axis=2,font.lab=2)
lines(Tobs$Control_Transpiration,Tobs$Control_Transpiration,lty=2,lwd=2)


qqplot(GWdata[,2],Out$gwlevel)
lines(c(-10,0),c(-10,0),lty=2,lwd=2)

# summarise the climate data
Climate <- data.frame(Dates = as.Date(Rain[,1]), Rain=Rain[,2],
	ETp = ETp[,2])

annual.clim <- aggregate(Climate[,2:3],list(year=substr(Climate$Dates,1,4)),sum)

barchart(ETp + Rain~year,annual.clim,col=c(1,"Gray50"),
		ylab=list("Total in mm", font=2,cex=1.3) ,
		scales=list(font=2,cex=1.3),
		key=list(text=list(lab=c("Potential ET","Rainfall"),cex=1.3),
        points=list(pch=15,col=c(1,"Gray50"),cex=1.3),
		x = .63, y = .92, corner = c(0, 0)))






##### older plotting ##########################

my.strip2 <- function(which.given, ..., factor.levels, bg) {
            levs <- c("Groundwater Transpiration",
				"Soil Transpiration/Evaporation", "Capillary flow")
            my.bg <- "gray90"
          strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}

bwplot(Tg*10  ~as.factor(loc),data=Out,subset=loc<100,
	xlab="Distance from river (m)",ylab="flux mm/day",outer=TRUE,
	strip=my.strip2,par.strip.text=list(font=2,cex=1.2))
xyplot(qcap*10~Dates|as.factor(loc),data=Out,subset=loc<100,
	xlab="Dates",ylab="flux mm/day",outer=TRUE,
	strip=my.strip,par.strip.text=list(font=2,cex=1.2),type="l")


bwplot(stress~as.factor(loc),data=Out,subset=loc<100,
	xlab="Distance from river (m)",ylab="flux mm/day",outer=TRUE,
	strip=my.strip2,par.strip.text=list(font=2,cex=1.2))

xyplot(qcap*10~Dates|as.factor(loc),data=Out,subset=loc<100,
	xlab="Dates",ylab="flux mm/day",outer=TRUE,
	strip=my.strip,par.strip.text=list(font=2,cex=1.2),type="l")



histogram(~(Ts*10)|as.factor(loc),data=Out,strip=my.strip,subset=loc<95,
	xlab="Soil transpiration (mm/day)",ylab="Frequency")
histogram(~(Ttotal*10)|as.factor(loc),data=Out,strip=my.strip,subset=loc<95,
	xlab="Total transpiration (mm/day)",ylab="Frequency")
histogram(~stress|as.factor(loc),data=Out,strip=my.strip, subset=loc<95,
	xlab="vegetation stress",ylab="Frequency")
#savePlot(paste(today,"HistTotalT_byLoc.jpeg",sep="_"),type="jpeg")
xyplot(stress~as.Date(Dates)|as.factor(loc),data=Out, subset=loc<100, 
	type="l",strip=my.strip,xlab = "Days", ylab = "static vegetation stress")
windows()
xyplot(Rain~as.Date(Dates)|as.factor(loc),data=Out, subset=loc<5, 
	type="l",strip=my.strip,xlab = "Days", ylab = "static vegetation stress")
bwplot(stress~as.factor(loc),data=Out,subset=loc<100,
	xlab="Distance from river (m)",ylab="dynamic stress",
	par.strip.text=list(font=2,cex=1.2))

xyplot(gwlevel~as.Date(Dates)|as.factor(loc),data=Out, subset=loc<100, 
	type="l",strip=my.strip,xlab = "Days", ylab = "groundwater height (m)")
histogram(~gwlevel|as.factor(loc),data=Out, subset=loc<100, 
 strip=my.strip,xlab = "groundwater height (m)")

qqplot(Out[Out$loc==0,"gwlevel"],Out[Out$loc==85,"gwlevel"])
lines(Out[Out$loc==0,"gwlevel"],Out[Out$loc==0,"gwlevel"],lty=2)


my.strip3 <- function(which.given, ..., factor.levels, bg) {
            levs <- paste("from river",c(0, 95, 525, 5100),"m")
            my.bg <- "gray90"
          strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}


sub_Out <- Out[(Out$loc %in% c(0, 95, 525,5100)),]

xyplot(gwlevel~as.Date(Dates)|as.factor(loc),data=sub_Out, 
	type="l",strip=my.strip3,xlab =list( "Days",font=2,cex=1.2), 
	ylab = list("groundwater height (m)",font=2,cex=1.2),
	par.strip.text=list(font=2,cex=1.2),ylim=c(-11,-3),
	scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),as.table=T)



#xyplot(gwlevel~as.Date(Dates)|as.factor(loc),data=Out, subset=loc<10, type="l",strip=my.strip,
#xlab = "Days", ylab = "groundwater height (m)",ylim=c(-6,-2))
finalgw <- result$gwlevel[1957,]
windows()
plot(distancetoriver,finalgw, type="l",
	xlab="Distance to river (m)",ylab="Groundwater depth below surface")#,
xlim=c(0,100))

save(finalgw,file="finalgwlevel.rdata")
plot(sort(Stream[,2]),ylim=c(0,0.1),type="l")

save.image(paste(today,"2dModelsession.rdata",sep="_"))
