# -------------------------------
# Simple example with uniform rainfall
# 

# to be removed
#-----------------------
require(lattice)
require(EcoHydro2D)
setwd("c:/users/rver4657/owncloud/ecohydro2d/examples")
#----------------------------

today <- format(Sys.Date(),"%Y%m%d")

# define uniform rainfall for two years
days <- seq(as.Date("2005-01-01"),as.Date("2006-12-31"),by=1)
Rain_uni <- data.frame(Dates = days, Rain = rep(0.1,730))
# sart with no water in Stream
Stream0 <- data.frame(Dates = days, Height = rep(0,730))
# and unifor ETp
ETp_uni <- data.frame(Dates = days, ET = rep(0.5,730))

# In advance, define all parameters in gwt_input_parameters.r
# check all input, see helpfile from GWT.exe for detail
edit(file="../input/gwt_input_parameters.r")

# Initial groundwater heads
gw_in <-  rep(-6,42)#(GWdata[1,2] + 0.5 - 0.048*(1:42))#seq(-8,-10,length=42) #finalgw 
# this is needed for Trees with deep roots
Zmean <- rep(600,42) 


# Define the vegetation series as being uniform grass
veg <- rep("Grass",42)
# Define what soil to use
soils <- "L Med Clay"
sp <- Soil_cpp("soils")
#
# No streamflow, just rainfall
# now includes separate K for aquifer (sand) based on Pfautsch et al.
result <- big.fun(N=nrow(Rain_uni),stype=soils,vtype=veg, 
                              aqK_in = sp$K_s/100,  aq_specy_in  = sp$spec_y,
                              Rain=Rain_uni, 
                  ETp=ETp_uni,stream=Stream0, gwheads = gw_in, Zmean = Zmean, 
                  today.m = today)

## Plotting example -------------------
# run stackfun on the output file
# you could also use reshape
Out <- stackfun(result,42)
head(Out) # creates a wide dataframe

# put in dates rather than sequential numbers:
Out$Dates <- rep(Rain_uni$Dates[1:length(result[[1]])],NX)

# Now use Lattice for some pretty plotting
# Set some background levels
trellis.par.set("regions", list(col=gray(25:100/100)))
trellis.par.set("background", list(col="white"))

#Locations to subset
unique(Out$loc)
Loc_choose <- c(0, 40, 80, 175)

# a fancy strip
my.strip3 <- function(which.given, ..., factor.levels, bg) {
            levs <- paste("from river",Loc_choose,"m")
            my.bg <- "gray90"
          strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}

# subset only a few locations (add 20 for river grid cell)
sub1_Out <- Out[(Out$loc %in% (Loc_choose+20)),]

xyplot(s + Ts*10~as.Date(Dates)|as.factor(loc),data=sub1_Out, 
	type="l",strip=my.strip3,xlab =list( "Days",font=2,cex=1.2), 
	ylab = list("Soil Saturation and Transpiration (cm)",font=2,cex=1.2),
	par.strip.text=list(font=2,cex=1.2), lty=c(1,2), col=c(1,"Gray50"),
	scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),as.table=T,
	key=list(text=list(lab=c("Soil","Transoiration"),cex=1),
        lines=list(lwd=2,lty=c(1,2),col=c(1,"Gray50")),x = .53, y = .75, corner = c(0, 0)))

# Now look at GWlevels
GWsim <- data.frame(Dates=sub1_Out$Dates, level=sub1_Out$gwlevel,
                    ind=rep("sim",nrow(sub1_Out)),loc=sub1_Out$loc)

# GW levels over time at selected locations
xyplot(level~as.Date(Dates)|as.factor(loc),data=GWsim,
	type="l",strip=my.strip3,xlab =list( "Days",font=2,cex=1.2), 
	ylab = list("height",font=2,cex=1.2), col=c(1:3,"gray50"),
	par.strip.text=list(font=2,cex=1.2),#ylim=c(-11,-3),
	scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),
	as.table=T,ylim=c(-10,5))

# also look at water table at a few dates
GW_dates <- as.Date(c("2005-07-01","2006-01-01","2006-07-01","2006-12-31"))
GWsim <- Out[(Out$Dates %in% GW_dates),]
GWsim <- GWsim[,c("Dates","gwlevel","gwrech","loc")]

my.strip3 <- function(which.given, ..., factor.levels, bg) {
  levs <- as.character(GW_dates)
  my.bg <- "gray90"
  strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}
# plot of river levels away from river at selected dates
xyplot(gwlevel~loc|as.factor(Dates),data=GWsim,
       type="l",strip=my.strip3,xlab =list( "m from river",font=2,cex=1.2), 
       ylab = list("height",font=2,cex=1.2), col=c(1:3,"gray50"),
       par.strip.text=list(font=2,cex=1.2),xlim=c(0,1000),
       scales=list(x = list(font=2,cex=1.2),y = list(font=2,cex=1.2)),
       as.table=T)

