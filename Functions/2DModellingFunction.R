# +++++++++++++++++++++++++++++++++++++++++++++++++ 
#  2 D modelling functions based on Joep vd Zanden's older work
#  rewrite to include the Rcpp calls
#  WV 20160628
# Can be sourced from a run scenario file
#----------------------------------------------
# General structure:
#
#	run gwt.exe to initialise
# Run loop through time * integration {
# 	Calculate all full days
#	if not a full day {
#		Calculate the water balance function (integration - 1) times
# 	else
#		Generate infiltration from rainfall function
#		if cumRecharge < threshold {
#			run water balance function (Rcpp)
#		else
#			run gwt.exe
#			run water balance function (Rcpp)
#		end if
#	end if
# end loop
#
#	run gwt.exe to complete
#
#	Produce output
#------------------------------------


# +++++++++++++++++++++++++++++++++++++++
# Function to run gwt.exe

Gwt.fun <- function(t,last_t_heads, BC,
				first=T, 
				GWcumvector,slopevector,
				aqK, aqy,
				DELT.m = DELT,
				NRBL.m = NRBL,
				river_in=list(length=4)) {

		# --------------------------------
				# Three steps in this function:
		# 1. write the gwt_timestepinput file using a function
		# 2. Execute GWT.exe
		# 3. Organise output and reset GWcumvector
		# -------------------------------------
    # variables in are:
    # t             : timestep value
    # last_t_heads  : groundwater heads at last time step
    # BC            : boundary conditions
    # first         : T or F, is model run for first time?
    # GWcumvector   : accumulated recharge since last run
    # slopevector   : Whether model has a slope and what the slope is
    # aqK           : Ks for aquifer (can be different from soil)
    # aqy           : specific yield for aquifer
    # DELT.m        : accumulated time since last run
    # NRBL.m/NRBL   : number of river blocks
    # river_in      : list of input with required river data (Ariver, criver, IXR and IXY)

		# 1. this writes the gwt_timestepinput file
	    # define input timestep
     	define_gwt_timestepinput(ITIM = 1, DELT = DELT.m, 
			  RECH = GWcumvector / DELT.m,
			  hriver = BC)
		# print(BC)
      # define the hydro input
  		define_gwt_hydroinput(last_t_heads, bottomvector,
			  ksat = aqK, spec_y = aqy)  #K_s from cm/d to m/d
  		# at first timestep define the grid and the river
		if (first == T) { # only for first timestep
	  		define_gwt_gridinput(NX, NY, NRBL.m, widthxvector=DELX, widthyvector=DELY)
		
			Ariver <- river_in[[1]]
			criver <- river_in[[2]]
			IXR <- river_in[[3]]
			IXY <- river_in[[4]]		
#		print(Ariver)
			if(NRBL.m > 0) {define_gwt_riverinput(NRBL.m,Ariver,criver,IXR,IXY)}
		}

		# 2. execute model
     	system("GWT.exe")
		
		# 3. organise output and reset values
		# All output should be prefixed by env$ to make sure they are 
		# reset DELT
     	DELT.m <- 1
		
		# recalculate heads
		heads <- matrixtovector(read.table("gwt_headoutput",fill = TRUE))[1:(NX * NY)]
		GWdepths <- heads - slopevector
    		delta_h <- heads - last_t_heads
     	qrivlat <- (delta_h - GWcumvector / aqy) * 100 * aqy
     	GWcum <- GWcumvector

       	## Include stop code: 
       	if(max(heads)> -1) {
         		finalN <- t
         		print(paste("GW reaches root zone for time step = ",finalN))
       # 		print(NameRun)
				print(GWcum)
				print(heads)
				print(DELT.m)
         		STOP<-TRUE
       	}

	# merge output
	out <- data.frame(heads = heads, GWdepths = GWdepths, delta_h = delta_h,
			qrivlat = qrivlat, Gwcum = GWcum)
	#print(out)
	return(list(DELT = DELT.m, out = out))
}
# ----------------------REVIEW BELOW for script
# To run and extract output:
# GWresult <- Gwt.fun(last_t_heads, first = F)
# DELT <- GWresult$DELT
# Results <- GWresult$out
# also reset GWcumvector
#  GWcumvector[1:(NX*NY)] <- 0 #Reset GWcumvector
# ------------------------

# +++++++++++++ end gw function ++++++++++++



# +++++++++++++++++++++++++++++++
# # Rainfall infiltration function
# Rain_in <- function(Ra ,n,Zr) { Ra/(n*Zr)}
# 
# # Ra is the rainfall vector
# this needs definition of "full day"
# +++++++++ end rainfall infiltration function

source("UtilityFunctions.R")



# +++++++++++++++++++++++++++++++++++++++++
# we are running over i timesteps and over j landscape cells
# GWthreshold = 0.0 m (set in gwt_input_parameters.txt)

big.fun <- function(N, stype, vtype, aqK_in, aq_specy_in, 
                    Rain, ETp, stream, gwheads, Zmean,
		rdir_in = "X:/vervoort/research/rcode/ecohydrology/2dmodelling",
		today.m = today, 
		fs=1, Res_in=100) {
#print(gwheads)

# -----------------------------
# Preliminaries
# N: number of days for run
# stype: soiltype from soilfunction.r
# vtype: vector of vegetation types for vegfunction.r length no of cells
# aqK_in and aq_specy_in s Aquifer K and sepcific yield
# Rain: 2 column dataframe: Date and Rainfall (cm/d)
# ETp: 2 column dataframe: Date and Pot ET (cm/d)
# adjusted WV 20120830
# ++++ inserted WV 20120912+++
# stream 3 column dataframe, second column is stream height, third column is right hand boundary
# +++++
##---Inserted WV 20121010 ---####
# gwheads is a scalar or vector of initial gwheads (negative in m) length = NX
# Zmean is a scalar of vector (length NX) of optimal groundwater depths for vegetation
# Zmean is in cm and positive (Can be -100*gwheads)
# rdir_in is the R directory to read functions and scripts
# --------------------------


  soilpar_in <- Soil_cpp(stype)
  n <- soilpar_in$n
  spec_y <- soilpar_in$spec_y

# read the stream and gw data  
  read.fun1(rdir = rdir_in,stream.m = stream, gwheads.m = gwheads, 
	  	Res = Res_in) #1/(soilpar_in$K_s*0.01))

	# create a name for any output files
	NameRun<-paste(today.m,N,stype,sep="_")

  # Initialise
  R <- Rain[,2] #vector with rain data
  deltat <- 1/12  
  DELT <- 1  
  finalN <- N #this will be overwritten when formula stops before N is reached  
  STOP <- FALSE
    # introduce slopevector
	if(NY==1){  slopevector.m <- slope
    			}else{
    		slopevector.m<-matrixtovector(slope)  }

# -------- end preliminaries ----------------------

# -------Create storage frames
	Storage_day <- list(id = 1:N, s = matrix(0,nrow=N,ncol=NX*NY), 
		qcap = matrix(0,nrow=N,ncol=NX*NY), Tgw = matrix(0,nrow=N,ncol=NX*NY),
		Tsoil = matrix(0,nrow=N,ncol=NX*NY), Ttotal = matrix(0,nrow=N,ncol=NX*NY), 
		Leakage = matrix(0,nrow=N,ncol=NX*NY), GWrech = matrix(0,nrow=N,ncol=NX*NY),
		smloss = matrix(0,nrow=N,ncol=NX*NY),static_stress = matrix(0,nrow=N,ncol=NX*NY), 
		surfoff = matrix(0,nrow=N,ncol=NX*NY))
# 	Storage_Gridday <- data.frame(s = numeric(length=NX*NY), 
# 		qcap = numeric(length=NX*NY), Tgw = numeric(length=NX*NY),
# 		Tsoil = numeric(length=NX*NY), Ttotal = numeric(length=NX*NY), 
# 		Leakage = numeric(length=NX*NY), GWrech = numeric(length=NX*NY),
# 		smloss = numeric(length=NX*NY),static_stress = numeric(length=NX*NY), 
# 		surfoff = numeric(length=NX*NY))
# 
# 	Storage_Subday <- data.frame(s = numeric(length=1/deltat), 
# 		qcap = numeric(length=1/deltat), Tgw = numeric(length=1/deltat),
# 		Tsoil = numeric(length=1/deltat), Ttotal = numeric(length=1/deltat), 
# 		Leakage = numeric(length=1/deltat), GWrech = numeric(length=1/deltat),
# 		smloss = numeric(length=1/deltat),static_stress = numeric(length=1/deltat), 
# 		surfoff = numeric(length=1/deltat))


	GW_store <- list(id = 1:(N+1), heads = matrix(0,nrow=N+1,ncol=NX*NY), 
				GWdepths = matrix(0,nrow=N+1,ncol=NX*NY),
			 	delta_h = matrix(0,nrow=N+1,ncol=NX*NY),	
				qrivlat = matrix(0,nrow=N+1,ncol=NX*NY), 
				GWcum = matrix(0,nrow=N+1,ncol=NX*NY))
# ----------------------------------------------
while (STOP == FALSE) {
# START OF THE LOOPS --------------------------#
# run a loop through the days
	for (t in 1:N) {
			# define initial soil moisture
		#	print(s_init) 
		if (t == 1)	s_init <- rep(0.5,NX*NY)
	  	# check the river level
		if (stream[t,2] < 0.06) {
			if (stream[t-1,2] >= 0.06 || t == 1) {
				#print("t = 1 or river dry")
			  # initialise gwheads
				if (t == 1) {last_heads_in <- gwheads} else {last_heads_in <- GW_store$heads[t-1,]} 
#				print(last_heads_in)
			  # initialise recharge
				if (t == 1) rech <- rep(0.0,NX) else rech <- GW_store$GWcum[t-1,]
               #  print(paste("Recharge=",rech))
				# run gwt.exe function first time to
				GW.out <- Gwt.fun(t, last_t_heads = last_heads_in, 
					BC = gwheads[NX], first=T, 
					aqK = aqK_in, aqy = aq_specy_in, GWcumvector = rech, 
					slopevector = slopevector.m,
					DELT.m = DELT,
					NRBL.m = 1,
					river_in = list(Ariver[2],criver[2],
					NX,NY))  
				# reset DELT (timestep)
				DELT <- GW.out$DELT
			} else { #this means river is dry, but last step also dry
				#print("t1 > 1")
				#print(GW_store$GWcum[t-1,])
				GW.out <- Gwt.fun(t,last_t_heads = GW_store$heads[t-1,], 
					first=F, BC = gwheads[NX], soilpar.m=soilpar_in, 
					aqK = aqK_in, aqy = aq_specy_in, GWcumvector=GW_store$GWcum[t-1,], 
					slopevector = slopevector.m,
					DELT.m = DELT,
					NRBL.m = 1)
			} #close loop whether first time run
		}  else  { # river level > 0
			# check if first time after river level back up or first time
			if (stream[t-1,2] < 0.06 || t == 1) {
				# run gwt.exe function
				if (t == 1) {
					last_heads_in <- gwheads
				} else {
					last_heads_in <- GW_store$heads[t-1,]
				}
				#print("river full") 
				#print(last_heads_in)
				if (t == 1) rech <- rep(0.0,NX) else rech <- GW_store$GWcum[t-1,]
		GW.out <- Gwt.fun(t,last_t_heads = last_heads_in, 
					first=T, BC = c(riverheads[t],gwheads[NX]),
					aqK = aqK_in,  aqy = aq_specy_in, GWcumvector=rech, 
					slopevector = slopevector.m,
					DELT.m = DELT,
					NRBL.m = NRBL,
					river_in = list(Ariver,criver,
					c(1,NX),c(1,NY)))  
			} else {
				GW.out <- Gwt.fun(t,last_t_heads = GW_store$heads[t-1,], 
   					first=F, BC = c(riverheads[t],gwheads[NX]),
					aqK = aqK_in,  aqy = aq_specy_in, GWcumvector=GW_store$GWcum[t-1,], 
					slopevector = slopevector.m,
					DELT.m = DELT,
					NRBL.m = NRBL)  
			} # Close if loop whether first time run
		} # close riverlevel loop
				# reset DELT
				#print(paste("gwt.exe is run and DELT =",DELT))
			#	DELT <- GW.out$DELT
				GW.out$out$Gwcum <- rep(0,NX)
			#	print(GW.out$out$Gwcum)
			#	print(GW_store$GWcum[t-1,])
			# change DELT
			DELT <- ifelse(GW.out$DELT==1,1,DELT + 1)
		# write away gw output 
		# if Gwt.fun wasn't run, it just writes the last one again
		r <- 1:(ncol(GW.out$out))
#		print(GW.out$out)
		GW_store <- list.write.fun(r, input=GW.out$out, 
			output = GW_store, t = t) 
		#str(GW_store)
#		print("GW_store$GWcum:")
#		print(GW_store$GWcum[t,])
#		save(GW_store,file="GW_store_temp")		
    Storage_Gridday <- WBEcoHyd(t=t,R=R,ET_in=ETp[,2],vtype=vtype,
                             soilpar=soilpar_in, s_init = s_init, 
                             fullday = t, Zmean = Zmean, 
                             GWdepths, GWdepths_prev, deltat=12, NX, NY)
		
		
 		# daily values across grid
		r <- 1:ncol(Storage_Gridday)
		Storage_day <- list.write.fun(r, input=Storage_Gridday, 
		 		output = Storage_day, t=t)
		#str(Storage_day)
		# Accumulate the GW recharge, this is now in GW_store
		GW_store$GWcum[t,] <- GW_store$GWcum[t,] + Storage_day$GWrech[t,]*0.01 # m/day		
	#	print(GW_store$GWcum[t,])
		print(paste("t =",t))
	} # close t loop

# ------------ END OF LOOPS ---------------------

# Run gwt_fun one last time
	GW.out <- Gwt.fun(t = N+1, last_t_heads = GW_store$heads[N,], first=F,#
				BC = c(riverheads[t],gwheads[NX]), 
				aqK = aqK_in,  aqy = aq_specy_in, 
				GWcumvector=GW_store$GWcum[N,], 
				slopevector = slopevector.m,
				DELT.m = DELT,
				NRBL.m = NRBL)
	# write away gw output 
	r <- 1:(ncol(GW.out$out))
	GW_store <- list.write.fun(r, input=GW.out$out, 
	 		output = GW_store, t = N + 1) 
STOP <- TRUE
}


# ---------------------------------------------
# Write and store output
# -------------------------------------
# write away values to some files (WB etc)
Storage.out <- do.call(cbind,Storage_day)
write.csv(Storage.out,paste("output/",NameRun,"_DailyOutput_table.csv",sep=""),
		row.names=FALSE)
GWstore.out <- do.call(cbind,GW_store)
write.csv(GWstore.out,paste("output/",NameRun,"_GWOutput_table.csv",sep=""),
		row.names=FALSE)

# Create a water balance table
# Not that important
Wb <- apply(Storage.out,2,sum)
Wb[c(1:42,((8*42):((9*42)-1)))] <- 
	apply(Storage.out,2,mean)[c(1:42,((8*42):((9*42)-1)))]
#print(Wb)

write.table(Wb,paste("output/",NameRun,"_Output_WB_table.csv",sep=""),row.names=FALSE,
    col.names=TRUE,sep=",")

# Create output out of function
# vector of X distances etc
Storage_day$x=distancetoriver
Storage_day$vtype = vtype
Storage_day$gwlevel = GW_store$GWdepths[1:N,]
Storage_day$qrivlat = GW_store$qrivlat[1:N,]
Storage_day$Rain = Rain[1:N,2]
Storage_day$Stream = Stream[1:N,2]
print(NameRun)
# ------------------------------------- end output-----------------
return(Storage_day)


} # end function
