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
				gwcumvector,slopevector_f,
				aq_K_f, aq_sy,
				DELT_f = DELT,
				river_f=list(length=4),
				gwinput_f=gwinput) {

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
    # gwcumvector   : accumulated recharge since last run
    # slopevector   : Whether model has a slope and what the slope is
    # aq_K_f           : Ks for aquifer (can be different from soil)
    # aq_sy           : specific yield for aquifer
    # DELT_f        : accumulated time since last run
    # river_f      : list of input with required river data (Ariver, criver, IXR and IXY)
    # gwinput_f     : list of input of required gwdata
      STOP = FALSE
    # define based on options
      NRBL = ecohydro2d.options()$NRBL
      NX = ecohydro2d.options()$NX
      NY = ecohydro2d.options()$NY
      DELX = ecohydro2d.options()$DELX
      DELY = ecohydro2d.options()$DELY
      
      # 1. this writes the gwt_timestepinput file
	    # define input timestep
     	define_gwt_timestepinput(ecohydro2d.options()$ITIM, DELT = DELT_f, 
			  RECH = gwcumvector / DELT_f,
			  hriver = BC)
		# print(BC)
      # define the hydro input
  		define_gwt_hydroinput(NX,NY,last_t_heads, gwinput_f$bottomvector,
			  ksat = aq_K_f, spec_y = aq_sy)  #K_s from cm/d to m/d
  		# at first timestep define the grid and the river
		if (first == T) { # only for first timestep
	  		define_gwt_gridinput(NX, NY, NRBL, 
	  		                     widthxvector=DELX, widthyvector=DELY)
		
			Ariver2 <- river_f[[1]]
			criver2 <- river_f[[2]]
			IXR2 <- river_f[[3]]
			IYR2 <- river_f[[4]]		
#		print(Ariver)
			if(NRBL > 0) {define_gwt_riverinput(NRBL,Ariver2,
			                                      criver2,IXR2,IYR2)}
		}

		# 2. execute model
  	  wd <- getwd()
  	  setwd(paste(wd,"gwtmod",sep="/"))
     	system("GWT.exe")
		  setwd(wd)
		# 3. organise output and reset values
		# All output should be prefixed by env$ to make sure they are 
		# reset DELT
     	DELT_f <- 1
		
		# recalculate heads
		heads <- matrixtovector(read.table("gwtmod/gwt_headoutput",
		                                   fill = TRUE))[1:(NX * NY)]
		gwdepths <- heads - slopevector_f
    		delta_h <- heads - last_t_heads
     	qrivlat <- (delta_h - gwcumvector / aq_sy) * 100 * aq_sy
     	gwcum <- gwcumvector

       	## Include stop code: 
       	if(max(heads)> -1) {
         		finalN <- t
         		print(paste("GW reaches root zone for time step = ",finalN))
       # 		print(NameRun)
# 				print(gwcum)
# 				print(heads)
# 				print(DELT_f)
         		STOP<-TRUE
       	}

	# merge output
	out <- data.frame(heads = heads, gwdepths = gwdepths, 
	                  delta_h = delta_h,
			qrivlat = qrivlat, gwcum = gwcum)
	#print(out)
	return(list(DELT = DELT_f, out = out, STOP = STOP))
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




# +++++++++++++++++++++++++++++++++++++++++
# we are running over i timesteps and over j landscape cells


big_fun <- function(N, stype, vtype, aq_K, aq_specy, 
                    Rain, ETp, stream, gwheads, Zmean=NULL,
		                today=format(Sys.Date(),"%Y%m%d")) {
#print(gwheads)
#browser()
# -----------------------------
# Preliminaries
# N: number of days for run
# stype: soiltype from soilfunction.r
# vtype: vector of vegetation types for vegfunction.r length no of cells
# aq_K_f and aq_specy_f is Aquifer K and sepcific yield
# Rain: 2 column dataframe: Date and Rainfall (cm/d)
# ETp: 2 column dataframe: Date and Pot ET (cm/d)
# stream 3 column dataframe, second column is stream height, third column is right hand boundary
# gwheads is a scalar or vector of initial gwheads (negative in m) length = NX
# Zmean is a scalar of vector (length NX) of optimal groundwater depths for vegetation
# Zmean is in cm and positive (Can be -100*gwheads)
# --------------------------


  soilpar_in <- Soil_cpp(stype)
  n <- soilpar_in$n
  # check if this does not interfer with aq specy, do we need this?
  spec_y <- soilpar_in$spec_y

# read the stream and gw data  
  if (!dir.exists(file.path("gwtmod"))) {
    print("GWT.exe and input files should be in a dir called GWTmod")
  } else {
    gwinput = read.fun1(stream.m = stream, gwheads.m = gwheads)
  }
# output of this function is a list with the following elements:
#   list(bottommatrix=bottommatrix,bottomvector=bottomvector,
#        slope=slope, riverheads=riverheads,
#        distancetoriver = distancetoriver,
#        init_heads = init_heads, criver = criver,
#        Ariver = Ariver, hriver=hriver)
   
  # definition of values from options()
  NX <- ecohydro2d.options()$NX
  NY <- ecohydro2d.options()$NY
  NRBL <- ecohydro2d.options()$NRBL
  IXR <- ecohydro2d.options()$IXR
  IYR <- ecohydro2d.options()$IYR
  
	# create a name for any output files
	NameRun<-paste(today,N,stype,sep="_")
  # Initialise
  R <- Rain[,2] #vector with rain data
  deltat <- 12  
  DELT <- 1  
  finalN <- N #this will be overwritten when formula stops before N is reached  
  STOP <- FALSE
    # introduce slopevector
	if(NY==1){  slopevector <- gwinput$slope
    			}else{
    		slopevector<-matrixtovector(gwinput$slope)  }

# -------- end preliminaries ----------------------

# -------Create storage frames
  
	Storage_day <- list(id = 1:N, s = matrix(0,nrow=N,ncol=NX*NY), 
		qcap = matrix(0,nrow=N,ncol=NX*NY), Tgw = matrix(0,nrow=N,ncol=NX*NY),
		Tsoil = matrix(0,nrow=N,ncol=NX*NY), Ttotal = matrix(0,nrow=N,ncol=NX*NY), 
		leakage = matrix(0,nrow=N,ncol=NX*NY), gwrech = matrix(0,nrow=N,ncol=NX*NY),
		smloss = matrix(0,nrow=N,ncol=NX*NY),static_stress = matrix(0,nrow=N,ncol=NX*NY), 
		surfoff = matrix(0,nrow=N,ncol=NX*NY))

# Storage for groundwater elements
	GW_store <- list(id = 1:(N+1), heads = matrix(0.0,nrow=N+1,ncol=NX*NY), 
				GWdepths = matrix(0.0,nrow=N+1,ncol=NX*NY),
			 	delta_h = matrix(0.0,nrow=N+1,ncol=NX*NY),	
				qrivlat = matrix(0.0,nrow=N+1,ncol=NX*NY), 
				GWcum = matrix(0.0,nrow=N+1,ncol=NX*NY))
# ----------------------------------------------
while (STOP == FALSE) {
# START OF THE LOOPS --------------------------#
pb <- txtProgressBar(min = 0, max = N, style = 3)
  
# run a loop through the days
	for (t in 1:N) {
#browser()
  	  # Initialise
  		if (t == 1)	{
  		  # 1. define initial soil moisture
  		  s_init <- rep(0.5,NX*NY)
    		# 2. initialise gwheads
  	   last_heads_in <- gwinput$init_heads
  	    # 3. initialise recharge
  	   rech <- rep(0.0,NX)
  		} else {
  	    last_heads_in <- GW_store$heads[t-1,]
  	    rech <- GW_store$GWcum[t-1,]
  	    #print(rech)
  		} 
  	  #browser()
  		# check the river level
  		if (stream[t,2] < 0.06) {
  			if (stream[t-1,2] >= 0.06 || t == 1) {
  				#print("t = 1 or river dry")
  				# run gwt.exe function first time to
  			  #browser()
  				GW.out <- Gwt.fun(t, last_t_heads = last_heads_in, 
  					BC = rep(gwinput$init_heads[length(gwinput$init_heads)],NRBL), first=T, 
  					aq_K_f = aq_K, aq_sy = aq_specy, gwcumvector = rech, 
  					slopevector_f = slopevector,
  					DELT_f = DELT,
  					river_f = list(gwinput$Ariver,gwinput$criver,IXR,IYR),
  					gwinput_f=gwinput)  
  				# reset DELT (timestep)
  				DELT <- GW.out$DELT
  			} else { #this means river is dry, but last step also dry
  				#print("t1 > 1")
  				#print(GW_store$GWcum[t-1,])
  				GW.out <- Gwt.fun(t,last_t_heads = GW_store$heads[t-1,], 
  					first=F, 
  					BC = rep(gwinput$init_heads[length(gwinput$init_heads)],NRBL),  
  					aq_K_f = aq_K, aq_sy = aq_specy, gwcumvector=GW_store$GWcum[t-1,], 
  					slopevector_f = slopevector,
  					DELT_f = DELT,
            gwinput_f=gwinput)
  			} #close loop whether first time run
  		}  else  { # river level > 0
  			# check if first time after river level back up or first time
  			if (stream[t-1,2] < 0.06 || t == 1) {
  				# run gwt.exe function
  			  #browser()
      		GW.out <- Gwt.fun(t,last_t_heads = last_heads_in, 
  					first=T, 
  					BC = c(gwinput$riverheads[t],gwinput$init_heads[length(gwinput$init_heads)]),
  					aq_K_f = aq_K,  aq_sy = aq_specy, gwcumvector=rech, 
  					slopevector_f = slopevector,
  					DELT_f = DELT,
  					river_f = list(gwinput$Ariver,gwinput$criver,IXR,IYR), 
  					gwinput_f=gwinput)
  			} else {
  				GW.out <- Gwt.fun(t,last_t_heads = GW_store$heads[t-1,], 
     					first=F, 
     					BC = c(gwinput$riverheads[t],gwinput$init_heads[length(gwinput$init_heads)]),
  					aq_K_f = aq_K,  aq_sy = aq_specy, gwcumvector=GW_store$GWcum[t-1,], 
  					slopevector_f = slopevector,
  					DELT_f = DELT,
  					gwinput_f=gwinput)  
  			} # Close if loop whether first time run
  		} # close riverlevel loop
  				# reset DELT
  				# print(paste("gwt.exe is run and DELT =",DELT))
  			#	DELT <- GW.out$DELT
  				GW.out$out$gwcum <- rep(0.0,NX)
  				STOP <- GW.out$STOP
  				#print(STOP)
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
  	#browser()
  		if(is.null(Zmean)==T) Zmean <- rep(200,NX) else Zmean <- Zmean 
      Storage_Gridday <- WBEcoHyd(t = t, R = R[t], ET_in = ETp[t,2],
                                vtype = vtype, soilpar = soilpar_in, 
                                s_init = s_init, Zmean = Zmean,
                                GWdepths = -100*GW_store$GWdepths[t,], 
                                GWdepths_prev = -100*last_heads_in,
                                deltat = deltat, NX, NY)
  		#browser()
   		# daily values across grid
  		Storage_day <- list.write.fun2(input = Storage_Gridday, 
  		 		output = Storage_day, t = t)
  		#str(Storage_day)
  		# Accumulate the GW recharge, this is now in GW_store
  		GW_store$GWcum[t,] <- GW_store$GWcum[t,] + Storage_day$gwrech[t,]*0.01 # m/day		
  		#print(Storage_day$GWrech[t,]*0.01)
  		  # update progress bar
  		  setTxtProgressBar(pb, t)
	} # close t loop
close(pb)
# ------------ END OF LOOPS ---------------------
#browser()
# Run gwt_fun one last time
	GW.out <- Gwt.fun(t = N+1, last_t_heads = GW_store$heads[N,], 
	      first=F,#
				BC = c(gwinput$riverheads[t],gwinput$init_heads[length(gwinput$init_heads)]), 
				aq_K_f = aq_K,  aq_sy = aq_specy, 
				gwcumvector=GW_store$GWcum[N,], 
				slopevector_f = slopevector,
				DELT_f = DELT,
				gwinput_f=gwinput)
	# write away gw output 
	r <- 1:(ncol(GW.out$out))
	GW_store <- list.write.fun(r, input=GW.out$out, 
	 		output = GW_store, t = N + 1) 
STOP <- TRUE
}


# ---------------------------------------------
# Write and store output
# -------------------------------------
# check if output dir exists
ifelse(!dir.exists(file.path(getwd(), "Output")), 
       dir.create(file.path(getwd(), "Output")), FALSE)
# write away values to some files (WB etc)
Storage.out <- do.call(cbind,Storage_day)
write.csv(Storage.out,paste(getwd(),"/output/",NameRun,"_DailyOutput_table.csv",sep=""),
		row.names=FALSE)
GWstore.out <- do.call(cbind,GW_store)
write.csv(GWstore.out,paste(getwd(),"/output/",NameRun,"_GWOutput_table.csv",sep=""),
		row.names=FALSE)

# Create a water balance table
# Not that important
Wb <- apply(Storage.out,2,sum)
Wb[c(1:42,((8*42):((9*42)-1)))] <- 
	apply(Storage.out,2,mean)[c(1:42,((8*42):((9*42)-1)))]
#print(Wb)

write.table(Wb,paste(getwd(),"/output/",NameRun,"_Output_WB_table.csv",sep=""),row.names=FALSE,
    col.names=TRUE,sep=",")

# Create output out of function
# vector of X distances etc
Storage_day$x=gwinput$distancetoriver
Storage_day$vtype = vtype
Storage_day$gwlevel = GW_store$GWdepths[1:N,]
Storage_day$qrivlat = GW_store$qrivlat[1:N,]
Storage_day$rain = Rain[1:N,2]
Storage_day$stream = stream[1:N,2]
#print(NameRun)
# ------------------------------------- end output-----------------
return(Storage_day)


} # end function
