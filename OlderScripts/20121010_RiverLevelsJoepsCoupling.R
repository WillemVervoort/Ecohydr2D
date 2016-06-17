today <- format(Sys.Date(),"%Y%m%d")
#setwd("C:/Documents and Settings/Administrator/Desktop/Hons/3-My Work/4 - R files")
setwd("x:/research/ecohydrology/2dmodelling")
#rdir <- "C:/Documents and Settings/Administrator/Desktop/Hons/3-My Work/4 - R files"
rdir <- "x:/research/rcode/ecohydrology/2dmodelling"

# In advance, define all parameters in gwt_input_parameters.txt

#### Read in the Rain data and climate data
# adjusted WV 20120830
Climate <- read.csv("YuleClim.csv")
Rain <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), Rain=Climate[,2]/10) # already in cm/d
ETp <- data.frame(dates=as.Date(Climate[,1], "%d/%m/%Y"), ETp=Climate[,5]/10) # translate to cm/day

# adjusted WV 20120928 
# monthly streamflow
Stream <- read.csv("20120928_yulerivlev.csv")

# Moving average on the stream flow
head(Stream)
Stream2 <- Stream
Stream2$Height <- filter(Stream$Height,filter=rep(0.02,50),method="convolution")

plot(as.Date(Stream$Date),Stream$Height,type="l")
lines(as.Date(Stream$Date),Stream2$Height,col="red")

# ------------------------------------------------------------------------------------
# ___________________________________________________________________________________
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Model: function of (N)
# ALL ONE BIG FUNCTION DOWN TO LINE 350

# Start function definition
Coupling <- function(N,stype,vtype,Rain,ETp,stream,gw_in,Zmean) {
	# N: number of days for run
	# stype: soiltype from soilfunction.r
	# vtype: vegetation type from vegfunction.r
	# Rain: 2 column dataframe: Date and Rainfall (cm/d)
	# ETp: 2 column dataframe: Date and Pot ET (cm/d)
	# adjusted WV 20120830
  # ++++ inserted WV 20120912+++
  # stream 2 column dataframe second column is stream height
  # +++++
  ##---Inserted WV 20121010 ---####
  # gw_in is a scalar or vector of initial gwheads (negative in m) length = NX
  # Zmean is a scalar of vector (length NX) of optimal groundwater depths for vegetation
  # Zmean is in cm and positive (Can be -100*gwheads)
  # --------------------------
  
	#read initial values and call other codes

	source(paste(rdir,"soilfunction.r",sep="/"))
	source(paste(rdir,"vegfunction.r",sep="/"))


	soilpar <- Soil(stype)
	vegpar <- Veg(vtype=vtype, soilpar=soilpar)
	DR <- vegpar$DR
  # This is a variable based on Vervoort and van der Zee (2012)
	vegpar$fs <- 0.2

	##### Adjusted by Willem 20121010
  # read in the gwheads and Zmean from input
  gwheads <-  gw_in
  # OLD c(-4,-5,-6,-6,-8,rep(-10.75,6)) # this is in m, but Zmean should be in cm and positive
	# OLD This is not very elegant, should actually be in vegpar, works for now
	# OLD Maybe Zmean should be higher than the reported "rest level" assuming trees closer to river might take up water better
	Zmean <- Zmean # OLD gwheads*-100


	# this attaching is needed to run the set of scripts below
	attach(list(gwheads=gwheads,stream=stream))

	source(paste(rdir,"gwt_input_parameters.r",sep="/")) #You need to adjust values in here as well!!
	source(paste(rdir,"20120724_FluxfunctionsforElise.R",sep="/"))  #!!!!! NEW means: adjusted formula.
	source(paste(rdir,"define_input.r",sep="/"))      #also the parameters are called

	# detach again to make sure there is no confusion
	detach(list(gwheads=gwheads,stream=stream))
	########### end adjustment Willem

	# create a name for any output files
	NameRun<-paste(today,N,stype,vtype,sep="_")


	# define initial river water level
	heads<-matrix(0,nrow=(N+1),ncol=(NX*NY)) 

	# initial heads
	heads[1,]<- gwheads#riverheads[1]
 
	if(length(DELX)>1){
  		est_leakage<-0.0*sum(Rain[,2])/nrow(Rain)  # start with analytical approximation. If 0, a flat gw table is used
  		total_xx <- 2*sum(DELX[2:NX])
  		for(ii in 2:length(DELX)){                 #est_leakage in cm/d. K_s in cm/d
    			heads[1,ii]<-sqrt(est_leakage/soilpar$K_s*(total_xx*distancetoriver[ii]-distancetoriver[ii]^2) +
      		(heads[1,ii]-bottom)^2)+bottom
    		}
  	#sum(DELX)
	}

  #first define initial values of BigFooFunction. Original vectors become matrices (1D to 2D)
  if (DR == TRUE) model <- "FB_new" else model <- rho_new_1 # "Fun_E" #
  sinit<-vegpar$s_star
  Zr <- vegpar$Zr
  n <- soilpar$n
  spec_y <- soilpar$spec_y
  R <- Rain[,2] #vector with rain data
  DELTvector<-R   # vector with DELT values to see when it is calculated again. Same length as Rain
  deltat <- 1/12    
  t<-1:(N*1/deltat)
  finalN <- N #this will be overwritten when formula stops before N is reached  
  ########
  D<-init_heads-bottom
 
  ss<-matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  Phi<-matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  Tg <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))  # Transpiration from deep roots
  qcap <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))  # capillary flux
  Ts <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY)) # Transpiration from soil
  Leakage <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))  # Leakage
  Tt <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))                                     
  smloss <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))                                     
  surfoff <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))                                     
  #LAI <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  #DM <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  #Aact <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  static_stress<-matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  ss[1,] <- sinit  # You need a begin value
  static_stress[1,]<- 0
  Tg[1,] <- 0
  Ts[1,] <- 0
  Leakage[1,]<-0
  smloss[1,] <- 0 
  #LAI[1,] <- vegpar$LAI
  #DM[1,] <- 0.5
  #Tt[1,] <- (1-exp(-0.45*LAI[1,]))*vegpar$Ep
  Tt[1,] <- (1-exp(-0.45*vegpar$LAI))*vegpar$Ep
  qcap[1,] <- 0
  #DORM <- FALSE
  DELTvector[1]<-1
  DELT <- 1
 
  # So let's now introduce the new parameters
  GWrech <- matrix(0,nrow=(N)*1/deltat,ncol=(NX*NY)) # Flux from groundwater (leakage - qcap - qzr)   
  GWrech[1,] <- 0
  GWcumvector <- vector("numeric",length=(NX*NY))
  GWcumvector[1:(NX*NY)]<-0
  GWcummatrix<-matrix(0,nrow=(N+1),ncol=(NX*NY)) 

  # introduce slopevector
  if(NY==1){  slopevector <- slope
    }else{
    slopevector<-matrixtovector(slope)  }

  #create vector with all h (=Z) values
  delta_h<-matrix(0,nrow=(N+1),ncol=(NX*NY)) 
  delta_h[1,] <- 0
  qrivlat<-matrix(0,nrow=(N+1),ncol=(NX*NY))   
  qrivlat[1,] <- 0
  GWdepths<-matrix(0,nrow=(N+1),ncol=(NX*NY)) 
  GWdepths[1,]<-gwheads-t(slopevector) # Might need to check whether this is right


#Define parameters in Toons model{
  #adjust parameters here if necessary: ksat, grid, recharge, ...
  #create grid. start with zero RECHarge vector for first run Toons model.
  define_gwt_gridinput(NX, NY, NRBL, widthxvector, widthyvector)#, NB, XB1, YB1, XB2, YB2, XB3, YB3, XB4, YB4)
  define_gwt_timestepinput(ITIM, DELT, RECH=GWrech[1,], hriver=riverheads[1])
  define_gwt_hydroinput(heads[1,], bottomvector,ksat=soilpar$K_s*0.01,spec_y=soilpar$spec_y)  #K_s from cm/d to m/d
  if(NRBL>0) {define_gwt_riverinput(rivermatrix=rivermatrix)}
 
 #Run Toons model once to get the proper head output / the first Z input
	#  system(paste("'",paste(rdir,"GWT.exe",sep="/"),"'"))
  	system("GWT.exe")

  #Loop: time This is the water balance with the GW model
  i<-1
  STOP<-FALSE  
  while(i<=(nrow(ss)-1)&STOP==FALSE) {
      # calculate whether rainfall infiltration occurs
     if (((floor(i*deltat)+1)==(i*deltat+1))&(max(abs(GWcumvector))<GWthreshold)&(DELT<DELTcrit)) {  #Means: every 12th time-step. 
      In <- R[floor(i*deltat)+1]/(n*Zr)#+ Flood[floor(i*deltat)+1]*deltat/(n*Zr)
       heads[floor(i*deltat)+1,]<-heads[floor(i*deltat),]
       GWdepths[floor(i*deltat)+1,]<-heads[floor(i*deltat)+1,]-t(slopevector)
       DELT<-DELT+1
     }  
     # evaluate cumulative GW flux, call GW model.
     # Calculate Toons model for the new time step. Not for each deltat time step, but for each time!!
     if (((floor(i*deltat)+1)==(i*deltat+1))&(max(abs(GWcumvector))>=GWthreshold)|(DELT==DELTcrit))  { 
	     In <- R[floor(i*deltat)+1]/(n*Zr)#+ Flood[floor(i*deltat)+1]*deltat/(n*Zr)

		## testing
		#print(paste(i,format(t(GWcumvector)/DELT,nsmall=5)))
		# this writes the gwt_timestepinput file
     	define_gwt_timestepinput(ITIM=(i*deltat+1), DELT=DELT, 
			RECH=format(t(GWcumvector)/DELT,nsmall=5),
			hriver=format(riverheads[floor(i*deltat)+1],nsmall=5))
		# execute model
     	system("GWT.exe")

     	DELT<-1
     	heads[floor(i*deltat)+1,]<-matrixtovector(read.table("gwt_headoutput",fill=TRUE))[1:(NX*NY)] #Fill: read also last lines when table is no square matrix
     	GWdepths[floor(i*deltat)+1,]<-heads[floor(i*deltat)+1,]-t(slopevector)
    		delta_h[floor(i*deltat)+1,]<- heads[floor(i*deltat)+1,]-heads[floor(i*deltat),]
     	qrivlat[floor(i*deltat)+1,]<-(delta_h[floor(i*deltat)+1,]-GWcumvector/spec_y)*100*spec_y
     	GWcummatrix[floor(i*deltat)+1,] <- t(GWcumvector)*100
     	GWcumvector[1:(NX*NY)]<-0 #Reset GWcumvector
       	## Include stop code: 
       	if(max(heads[floor(i*deltat)+1,])>-Zr/100) {
         		finalN <- i*deltat
         		print(paste("GW reaches root zone for time step = ",finalN))
         		print(NameRun)
         
         		#      heads[N+1,]<-(-Zr/100)
         		STOP<-TRUE
       	}
       } 
     if (((floor(i*deltat)+1)!=(i*deltat+1))) {  #Means: every time step except the 12th  
       In <- 0
     }
     # calculate root zone model for each Z in Toons model --> Loop: space
     if(STOP==FALSE){
     for(j in 1:(NX*NY)) {
        # define Z values 
        Z <- GWdepths[floor(i*deltat)+1,j]*-100
        # The Big Foo Loop pasted!
        #  Empty[i]<-n*Zr*(1-ss[i])

        Phi[i,j]<-In
 ## Testing
 #       print(Phi[i,j])
 #       print(ss[i,j])
        if(Phi[i,j]>(1-ss[i,j])){
          Phi[i,j]<- 1-ss[i,j]
          surfoff[i+1,j]<- (In-Phi[i,j])*n*Zr
          }
      # in here read in the actual pot ET into vegpar$Ep
        ####### inserted WV 20120830
        vegpar$Ep <- ETp[floor(i*deltat)+1,2]
        loss <- do.call(model,list(s=ss[i,j],soilpar=soilpar,vegpar=vegpar,
                        Z=Z,Zmean=Zmean[j],Z.prev=(GWdepths[max(1,floor(i*deltat)),j]*-100))) #i should be at least 1
        #########
        # update the water balance
# ss.t is temporary ss
      ss.t <- ss[i,j] + Phi[i,j] -loss[1]/(n*Zr)*deltat
      # build in safety if capillary flux fills soil with dynamic water table
      if (ss.t > 1) {  # very wet conditions
          ss[i+1,j] <- 1
#          qcap[i+1,j] <- (1 - ss[i,j])*n*Zr/deltat
          qcap[i+1,j] <- loss[4] - (ss.t - 1)*n*Zr/deltat
      } else { # normal conditions
          ss[i+1,j]<-ss.t
          qcap[i+1,j] <- loss[4]#/(n*Zr)*deltat
      }
        static_stress[i+1,j]<- zeta(s=ss[i+1,j],s_star=vegpar$s_star,s_w=vegpar$s_w,q=q)
        Tg[i+1,j] <- loss[3]#/(n*Zr)*deltat
 #       print(loss[3])
 #       qcap[i+1,j] <- loss[4]#/(n*Zr)*deltat           ###!!!
        Ts[i+1,j] <- loss[2]#/(n*Zr)*deltat
        Tt[i+1,j] <- Tg[i+1,j] + Ts[i+1,j]
        Leakage[i+1,j] <- ifelse(loss[5]>0,loss[5],0)

        smloss[i+1,j] <- loss[1]  - (qcap[i+1,j] - loss[4]) #increase with difference between originally calculated and new (pot.reduced) capillary flux
        GWrech[i+1,j] <- Leakage[i+1,j]-loss[3]-qcap[i+1,j]

# Accumulate the GW recharge
     GWcumvector[j] <-GWcumvector[j] + GWrech[i+1,j]*deltat*0.01 #so in m/d
     }
  # testing
  #print(i)
  #print(ss[i+1,j])
    
  #	print(heads[i,])
  i<-i+1
 }}
# execute GWT model one last time. + create final heads
   define_gwt_timestepinput(ITIM=(N+1), DELT=DELT, 
		RECH=format(t(GWcumvector)/DELT,nsmall=5),
		hriver=format(riverheads[floor(i*deltat)+1],nsmall=5))
   system("GWT.exe")     
   heads[N+1,]<-matrixtovector(read.table("gwt_headoutput",fill=TRUE))[1:(NX*NY)] #Fill: read also last lines when table is no square matrix
   GWdepths[N+1,]<-heads[N+1,]-t(slopevector) #heads & GWdepths N+1 are not in output table. Just used to calculate dh for waterbalance.
   GWcummatrix[N+1,]<-t(GWcumvector)*100
   delta_h[N+1,]<- heads[N+1,]-heads[N,]    
   qrivlat[N+1,]<-(delta_h[N+1,]-GWcumvector/spec_y)*100*spec_y

                       

  # create output data
        days <- sort(rep(1:N,1/deltat))
        sat_out <- aggregate(ss,list(day=days),mean)
        static_stress_out<-aggregate(static_stress,list(day=days),mean)
        Ts_out <- aggregate(Ts,list(day=days),mean) # Transpiration from soil bucket
        Tg_out <- aggregate(Tg,list(day=days),mean)# Transpiration from the groundwater
#        print(Tg_out)
        smloss_out <- aggregate(smloss,list(day=days),mean) # Actual loss of soil moisture
        surfoff_out <- aggregate(surfoff,list(day=days),mean) # surface runoff
        qcap_out <- aggregate(qcap,list(day=days),mean) # capillary flow
        leak_out <- aggregate(Leakage,list(day=days),mean) # Deep drainage
        GWrech_out <- aggregate(GWrech,list(day=days),mean)  # back to cm/d
#print(GWdepths)
# write important values to hard disk (s, Z, fluxes Q, ET)
#create output_table, start with line with initial values
  #write.table(initial_table,paste(rdir,"Output_table",sep="/"),row.names=FALSE,col.names=TRUE,sep=" ")
#create table with output results       
output_table<- data.frame(days=sat_out[,1], Rain=R[1:N], River=riverheads[1:N],
	Heads=heads[1:N,1:(NX*NY)], GWdepth=GWdepths[1:N,1:(NX*NY)],
     smloss=smloss_out[,2:(NX*NY+1)],Ts=Ts_out[,2:(NX*NY+1)],  Tg=Tg_out[,2:(NX*NY+1)], 
    qcap=qcap_out[,2:(NX*NY+1)], leak=leak_out[,2:(NX*NY+1)], GWrech=GWrech_out[,2:(NX*NY+1)],
    GWrechcum=GWcummatrix[2:(N+1),1:(NX*NY)], deltaH=delta_h[2:(N+1),1:(NX*NY)],surfoff=surfoff_out[,2:(NX*NY+1)], 
    qrivlat=qrivlat[2:(N+1),1:(NX*NY)], s=sat_out[,2:(NX*NY+1)], st_stress=static_stress_out[,2:(NX*NY+1)])
#  LAI=LAI.out[,2:(NX*NY+1)], DM=DM.out[,2:(NX*NY+1)], Aact=Aact.out[,2:(NX*NY+1)],
#paste output results in output_table
write.csv(output_table,paste(NameRun,"Output_table.csv",sep="_"),row.names=FALSE)
     

# Create a water balance table

#deltah<-(heads[N+1,]-heads[1,])
sum_rain<-sum(R[2:finalN])
sum_smloss<-vector(mode="numeric",length=NX*NY)
sum_Ts<-vector(mode="numeric",length=NX*NY)
sum_Tg<-vector(mode="numeric",length=NX*NY)
sum_qcap<-vector(mode="numeric",length=NX*NY)
sum_leak<-vector(mode="numeric",length=NX*NY)
sum_GWrech<-vector(mode="numeric",length=NX*NY)
sum_surfoff<-vector(mode="numeric",length=NX*NY)
sum_delta_h<-vector(mode="numeric",length=NX*NY)
deltas_calc<-vector(mode="numeric",length=NX*NY)
deltas_real<-vector(mode="numeric",length=NX*NY)
mean_s.T<-vector(mode="numeric",length=NX*NY)
mean_s.2<-vector(mode="numeric",length=NX*NY)
mean_Z.T<-vector(mode="numeric",length=NX*NY)
mean_Z.2<-vector(mode="numeric",length=NX*NY)
var_s.T<-vector(mode="numeric",length=NX*NY)
var_s.2<-vector(mode="numeric",length=NX*NY)
var_Z.T<-vector(mode="numeric",length=NX*NY)
var_Z.2<-vector(mode="numeric",length=NX*NY)
cell<-vector(mode="numeric",length=NX*NY)
deltah_GWrech<-vector(mode="numeric",length=NX*NY)
sum_qrivlat<-vector(mode="numeric",length=NX*NY) # river and lateral inflow per cell. Difference between deltah and deltah_GWrech

for(i in 1:(NX*NY)){
  sum_smloss[i]<-sum(smloss_out[2:finalN,i+1])
  sum_Ts[i]<-sum(Ts_out[2:finalN,i+1])
  sum_Tg[i]<-sum(Tg_out[2:finalN,i+1])
  sum_qcap[i]<-sum(qcap_out[2:finalN,i+1])
  sum_leak[i]<-sum(leak_out[2:finalN,i+1])
  sum_GWrech[i]<-sum(GWrech_out[1:finalN,i+1]) 
  sum_surfoff[i]<-sum(surfoff_out[2:finalN,i+1])
  sum_delta_h[i]<-sum(delta_h[2:(N+1),i])
  deltas_real[i]<- sat_out[finalN,i+1]-ss[1,i]
  deltas_calc[i]<-(sum_rain-sum_smloss[i]-sum_surfoff[i])/(vegpar$Zr*n)
  mean_s.T[i] <- mean(sat_out[2:finalN,i+1])            # calculate mean_s
  mean_s.2[i] <- mean(sat_out[(floor(N/2)):finalN,i+1]) # calculate mean_s of second half of run - check if the mean changes.
  mean_Z.T[i] <- mean(heads[1:(N+1),i])            # calculate mean_Z
  mean_Z.2[i] <- mean(heads[(floor(N/2)):(N+1),i]) # calculate mean_Z of second half of run - check if the mean changes.
  var_s.T[i] <- var(sat_out[2:finalN,i+1])            # calculate var_s
  var_s.2[i] <- var(sat_out[(floor(N/2)):finalN,i+1]) # calculate var_s of second half of run - check if the mean changes.
  var_Z.T[i] <- var(heads[1:(N+1),i])            # calculate var_Z
  var_Z.2[i] <- var(heads[(floor(N/2)):(N+1),i]) # calculate var_Z of second half of run - check if the mean changes.
  cell[i]<-i
  sum_qrivlat[i]<-sum(qrivlat[2:(N+1),i])             
 }

# Write and store output
WB_table<- data.frame(cell=cell, slope=slopevector, deltah=sum_delta_h, 
    sum_rain=sum_rain, sum_smloss=sum_smloss, sum_Ts=sum_Ts, sum_Tg=sum_Tg, sum_qcap=sum_qcap, 
    sum_leak=sum_leak, sum_GWrech=sum_GWrech, sum_qrivlat=sum_qrivlat, sum_surfoff=sum_surfoff, 
    deltas_real=deltas_real, deltas_calc=deltas_calc, mean_s.T=mean_s.T, mean_s.2=mean_s.2,
    var_s.T=var_s.T, var_s.2=var_s.2, mean_Z.T=mean_Z.T, mean_Z.2=mean_Z.2, var_Z.T=var_Z.T,
    var_Z.2=var_Z.2) 
write.table(WB_table,paste(NameRun,"Output_WB_table.csv",sep="_"),row.names=FALSE,
    col.names=TRUE,sep=" ")


print(NameRun)

#Output for stack formula:
return(data.frame(x=distancetoriver,heads0=heads[1,],headshw=heads[floor(N/2),],
    headsN=heads[N+1,],meanZ=mean_Z.T,varZ=var_Z.T,mean_s=mean_s.T,var_s=var_s.T,
    meanZ2=mean_Z.2,varZ2=var_Z.2,mean_s2=mean_s.2,var_s2=var_s.2,
    VEG=vtype,SOIL=stype,K_s=soilpar$K_s))
}
# END OF COUPLING FUNCTION
# ------------------------------------------------------------------------------------
# ___________________________________________________________________________________
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Now run a few different simulations
# Inserted 20121010 WV
# define initial gwheads and Zmean
gw_in <- c(-4,-5,-6,-6,-8,rep(-10.75,6))
Zmean <- rep(500,11)  # gw_in*-100 could also try Zmean <- 500


# TO MODEL DEEP ROOTS WITH TREES YOU NEED TO CHANGE vtype to
# vtype = "TreesDR" otherwise "TreesNoDR"
Coupling(N=250,stype="L Med Clay", vtype="TreesDR",Rain=Rain, 
         ETp=ETp,stream=Stream,gw_in = gw_in, Zmean = Zmean)

# assign to look at output
result <- Coupling(N=nrow(Rain),stype="L Med Clay", vtype="TreesDR",
                  Rain=Rain,ETp=ETp, stream=Stream,gw_in = gw_in, Zmean = Zmean)
result$x # should give you distance to river


##### SOme test plotting, example by Willem ##################

# look at some output
par(mfrow=c(1,1))
Outfile <- read.csv("20121011_1957_L Med Clay_TreesDR_Output_table.csv")
par(mar=c(4,4,2,5)) # Leave space for z axis
plot(Outfile$days,Outfile$Rain, type="h",
     xlab = "Time (days)", ylab="Rain (cm/d),Stage (m)")
lines(Outfile$days,(Outfile$River+3),col="purple",lwd=2)
par(new = TRUE)
plot(Outfile$days,Outfile$Heads.1, type = "l", axes=FALSE, bty = "n", 
	xlab = "", ylab="",col="red",ylim=c(-11,-3))
axis(side=4, at = pretty(range(Outfile[,grep("Heads",colnames(Outfile))])))
mtext("Groundwater depth (m)", side=4, line=3)
lines(Outfile$days,Outfile$Heads.2,lty=2,col="red")
lines(Outfile$days,Outfile$Heads.3,lty=3,col="red")
lines(Outfile$days,Outfile$Heads.4,lty=4,col="red")
lines(Outfile$days,Outfile$Heads.5,lty=1,col="blue")
lines(Outfile$days,Outfile$Heads.6,lty=2,col="blue")
lines(Outfile$days,Outfile$Heads.7,lty=3,col="blue")
lines(Outfile$days,Outfile$Heads.8,lty=4,col="blue")
lines(Outfile$days,Outfile$Heads.9,lty=1,col="green")
lines(Outfile$days,Outfile$Heads.10,lty=2,col="green")
lines(Outfile$days,Outfile$Heads.11,lty=3,col="green")

lgd.txt <- c("Rain","River stage",paste("heads",cumsum(DELX),"m"))
legend("topleft",lgd.txt,col=c(1,"purple",rep("red",4),rep("blue",4),rep("green",3)),
	lty=c(1,1,1:4,1:4,1:3),lwd=c(1,2,rep(1,11)),cex=0.7)

# take out the Ts and Tg values
Ts <- Outfile[,grep("Ts",colnames(Outfile))]
Tg <- Outfile[,grep("Tg",colnames(Outfile))]
TT <- Ts + Tg
par(mfrow=c(3,4))
for (i in 1:ncol(TT)) {
  qqplot(Climate[,3],TT[,i]*10,xlab="observed ET (mm/day)",ylab="simulated ET (mm/day)")
  abline(0,1)
  legend("topleft",paste("loc",i),pch=1)
}




## More pretty plotting? example by Willem ####-------------------
require(lattice)

# Do a stack operation, Joep actually wrote a function for this
names <- colnames(Outfile)
stackfun <- function(data,NX) { # this needs to be adapted relative to what you want
	df <- data.frame(days=rep(data$days,NX), Rain=rep(data$Rain,NX),
			River=rep(data$River,NX), Ts = stack(data[,grep("Ts",names)]),
			Tg = stack(data[,grep("Tg",names)]), 
			s = stack(data[,grep("s.V",names,fixed=T)[23:33]]),
			stress = stack(data[,grep("stress",names)]))
	return(df)
}
# run stackfun on the output file
Out <- stackfun(Outfile,11)
head(Out)

# Now use Lattice for some pretty plotting
#Define a strip function
my.strip <- function(which.given, ..., factor.levels, bg) {
            levs <- paste("From river",cumsum(c(20,30,50,100,100,200,200,300,500,1000,5000)),"m")
            my.bg <- "gray90"
          strip.default(which.given, ..., factor.levels = levs, bg = my.bg)
}

# Set some background levels
trellis.par.set("regions", list(col=gray(25:100/100)))
trellis.par.set("background", list(col="white"))
# Make a histogram of the s values
histogram(~s.values|s.ind,data=Out,strip=my.strip,
	xlab="Soil saturation",ylab="Frequency")

# Make a histogram of Total Transpiration to compare to observed
histogram(~(Ts.values+Tg.values)*10|Ts.ind,data=Out,strip=my.strip,
	xlab="Total transpiration (mm/day)",ylab="Frequency")
savePlot(paste(today,"HistTotalT_byLoc.jpeg",sep="_"),type="jpeg")


