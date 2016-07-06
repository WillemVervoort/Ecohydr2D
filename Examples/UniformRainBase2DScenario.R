# -------------------------------
# Some real scenario runs
# 
require(lattice)

require(EcoHydro2D)

#setwd("x:/vervoort/research/ecohydrology/2dmodelling")
setwd("c:/users/rver4657/owncloud/ecohydrology2dmodellng/examples")

today <- format(Sys.Date(),"%Y%m%d")

# read in all data
load("Allinput.rdata")

# define uniform rainfall
Rain_uni <- Rain
Rain_uni[,2] <- rep(0.1,nrow(Rain))
Stream0 <- Stream
Stream0[,2] <- 0
ETp_uni <- ETp
ETp_uni[,2] <- 0.5

# In advance, define all parameters in gwt_input_parameters.r
#edit(file="gwt_input_parameters.r")
# gives error

# Now run a few different simulations
# define initial gwheads and Zmean
# gw_in needs to be positive and in cm
gw_in <-  rep(-6,42)#(GWdata[1,2] + 0.5 - 0.048*(1:42))#seq(-8,-10,length=42) #finalgw 
Zmean <- rep(600,42) 
#fs_veg <- 0.25 # this is the spreading of the root water uptake around z_mean

# control the randonmness
#set.seed(1000)

# Define the vegetation series as being uniform grass
veg <- rep("Grass",42)
# Insert some bare soil in the series
#veg[sample(1:length(veg),size=ceiling(0.3*length(veg)))] <- 
#			"Bare"
# define the trees section
#veggies <- c("Bare",rep("TreesDR",10), veg,"Bare")
# Define what soil to use
# in the future this should become also a vector
soils <- "L Med Clay"
sp <- Soil_cpp("Loamy Sand")
# define spec_y separately
sp$spec_y <- 0.5
#
# No streamflow, just rainfall
# now includes separate K for aquifer (sand) based on Pfautsch et al.
system.time({result <- big.fun(N=500,stype=soils,vtype=veg, 
                              aqK_in = sp$K_s/100,  aq_specy_in  = sp$spec_y,
                              Rain=Rain_uni, 
                  ETp=ETp_uni,stream=Stream0, gwheads = gw_in, Zmean = Zmean, 
                  today.m = today, fs = fs_veg)})

## Plotting example by Willem ####-------------------


# run stackfun on the output file
Out <- stackfun(result,42)
head(Out) # creates a wide dataframe

# put in dates rather than sequential numbers:
Out$Dates <- rep(Rain$dates[1:length(result[[1]])],42)

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

my.strip3 <- function(which.given, ..., factor.levels, bg) {
            levs <- paste("from river",c(0, 40, 80, 175),"m")
            my.bg <- "gray90"
          strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}


sub1_Out <- Out[(Out$loc %in% c(20, 60, 100, 295)),]

xyplot(s+Ts*10~as.Date(Dates)|as.factor(loc),data=sub1_Out, 
	type="l",strip=my.strip3,xlab =list( "Days",font=2,cex=1.2), 
	ylab = list("Transpiration (cm)",font=2,cex=1.2),
	par.strip.text=list(font=2,cex=1.2), lty=c(1,2), col=c(1,"Gray50"),
	scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),as.table=T,
	key=list(text=list(lab=c("Groundwater","Soil"),cex=1),
        lines=list(lwd=2,lty=c(1,2),col=c(1,"Gray50")),x = .53, y = .75, corner = c(0, 0)))


GWsim <- data.frame(Dates=sub1_Out$Dates, level=sub1_Out$gwlevel,
                    ind=rep("sim",nrow(sub1_Out)),loc=sub1_Out$loc)
GWobs <- data.frame(Dates=rep(as.Date(GWdata$Date,"%d/%m/%Y"),4), 
                    level=rep(GWdata$Control_m_f_surf,4),
                    ind=rep("obs",nrow(GWdata)*4),
                    loc=rep(c(20, 60, 100, 295),each=nrow(GWdata)))
Strobs <- data.frame(Dates=rep(as.Date(Stream$Date,"%d/%m/%Y"),4), 
                    level=rep(Stream$Height,4),
                    ind=rep("stream",nrow(Stream)*4),
                    loc=rep(c(20, 60, 100, 295),each=nrow(Stream)))
#GWplot <- rbind(GWsim,GWobs,Strobs,Rainobs)
GWplot <- rbind(GWsim,GWobs,Strobs)


xyplot(level~as.Date(Dates)|as.factor(loc),data=GWplot, groups=ind,
	type="l",strip=my.strip3,xlab =list( "Days",font=2,cex=1.2), 
	ylab = list("height",font=2,cex=1.2), col=c(1:3,"gray50"),
	par.strip.text=list(font=2,cex=1.2),#ylim=c(-11,-3),
	scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),
	as.table=T,ylim=c(-10,5))

