# +++++++++++++++++++++++++++++++++++++++++++++++++ 
#  2 D modelling functions based on Joep vd Zanden's older work
#  WV 20130615
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
#			run water balance function
#		else
#			run gwt.exe
#			run water balance function
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

Gwt.fun <- function(t,last_t_heads, RB = gwheads[NX],first=T, 
				GWcumvector,slopevector,
				soilpar.m = soilpar,
				DELT.m = DELT,
				NRBL.m = NRBL) {

		# --------------------------------
		# This function needs a specified environment
		# Gwt.fun(Myenv)
		# Three steps in this function:
		# 1. write the gwt_timestepinput file using a function
		# 2. Execute GWT.exe
		# 3. Organise output and reset GWcumvector
		# -------------------------------------

		# 1. this writes the gwt_timestepinput file
	
     	define_gwt_timestepinput(ITIM = t, DELT = DELT.m, 
			RECH = GWcumvector / DELT.m,
			hriver = c(riverheads[i],RB))
		if (first == T) { # only for first timestep
	  		define_gwt_gridinput(NX, NY, NRBL.m, widthxvector=DELX, widthyvector=DELY)#, NB, XB1, YB1, XB2, YB2, XB3, YB3, XB4, YB4)
	  		define_gwt_hydroinput(last_t_heads, bottomvector,
				ksat = soilpar.m$K_s * 0.01, spec_y = soilpar.m$spec_y)  #K_s from cm/d to m/d
  	
			if(NRBL.m > 0) {define_gwt_riverinput(rivermatrix = rivermatrix)}
		}

		# 2. execute model
     	system("GWT.exe")
		
		# 3. organise output and reset values
		# All output should be prefixed by env$ to make sure they are 
		# reset DELT
     	DELT.m <- 0
		
		# recalculate heads
		heads <- matrixtovector(read.table("gwt_headoutput",fill = TRUE))[1:(NX * NY)]
		GWdepths <- heads - slopevector
    		delta_h <- heads - last_t_heads
     	qrivlat <- (delta_h - GWcumvector / soilpar.m$spec_y) * 100 * soilpar.m$spec_y
     	GWcum <- GWcumvector * 100

       	## Include stop code: 
       	if(max(heads)> -1) {
         		finalN <- t
         		print(paste("GW reaches root zone for time step = ",finalN))
         		print(NameRun)
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
# Rainfall infiltration function
Rain_in <- function(fullday,R ,n,Zr) { R[fullday]/(n*Zr)}

# R is the rainfall vector
# this needs definition of "full day"
# +++++++++ end rainfall infiltration function


# +++++++++++++++++++++++++++++++++++++++++++++
# Water balance function
# Just calculates the water balance function once
WB_fun <- function(vegpar.m, In, fullday, last_t_soil.sat, 
			soilpar.m=soilpar, 
			model = model_in, 
			Zmean.m,
			deltat.m = deltat,
		   	Z.prev = (GWdepths[max(1,(fullday-1)),j]*-100),
		   	Z.m = GWdepths[fullday,j]*-100) {
		

		n <- soilpar.m$n
		Zr <- vegpar.m$Zr

		# Rainfall into infiltration
        Phi<-In
		# check if soil bucket full, create runoff
        if(Phi > (1 - last_t_soil.sat)){
          Phi <- 1 - last_t_soil.sat
          surfoff <- (In - Phi)*n*Zr
          } else { surfoff <- 0 }
		# Call the loss model
         if (vegpar.m$DR==T) {
		loss <- do.call(model,list(s = last_t_soil.sat, 
					  soilpar = soilpar.m, vegpar = vegpar.m,
                        Z = Z.m, Zmean = Zmean.m, Z.prev = Z.prev))
		} else {
		loss <- do.call(model,list(s = last_t_soil.sat, Z = Z.m,
   					  soilpar = soilpar.m, vegpar = vegpar.m,
                        Z.prev = Z.prev))
		}
        # update the water balance
		# ss.t is temporary ss
     	 ss.t <- last_t_soil.sat + Phi - loss[1]/(n*Zr)*deltat.m
	      # build in safety if capillary flux fills soil with dynamic water table
	          static_stress<- zeta(s=ss.t,s_star=vegpar.m$s_star,s_w=vegpar.m$s_w,q=vegpar.m$q)
    if (ss.t > 1) {  # very wet conditions
    		     soil.sat <- 1
          	 qcap <- loss[4] - (ss.t - 1)*n*Zr/deltat.m
      	  } else { # normal conditions
          	 soil.sat<-ss.t
          	 qcap <- loss[4]
      	  }
 	       Tg <- loss[3]
 #       print(loss[3])
	        Ts <- loss[2]
    		    Tt <- Tg + Ts
        	    Leakage <- ifelse(loss[5]>0,loss[5],0)

	        smloss <- loss[1]  - (qcap - loss[4]) #increase with difference between originally calculated and new (pot.reduced) capillary flux
    		    GWrech <- Leakage-loss[3]-qcap

		out <- list(s = soil.sat, qcap = qcap, Tgw = Tg, Tsoil = Ts, Ttotal = Tt, 
				Leakage = Leakage, GWrech = GWrech, smloss = smloss, 
				static_stress = static_stress, surfoff = surfoff)

}

# Need to run this function and store output somewhere in a daily frame


# Read in function for sourcing other scripts
read.fun <- function(rdir = rdir_in, stream.m = stream, gwheads.m = gwheads) {
			source(paste(rdir,"soilfunction.r",sep="/"))
			source(paste(rdir,"vegfunction.r",sep="/"))
			attach(list(gwheads=gwheads.m,stream=stream.m))

			source(paste(rdir,"gwt_input_parameters.r",sep="/")) #You need to adjust values in here as well!!
			source(paste(rdir,"20120724_FluxfunctionsforElise.R",sep="/"))  #!!!!! NEW means: adjusted formula.
			source(paste(rdir,"define_input.r",sep="/"))      #also the parameters are called

			# detach again to make sure there is no confusion
			detach(list(gwheads=gwheads.m,stream=stream.m))
}

# Create a write function for storage matrices
list.write.fun <- function(r,output,input,t) {
	
	# input and output are a dataframe and a list
	# i is the time counter
	# r is a counter along the columns of inputs
	for (r in 1:length(r)) {
		output[[names(output)[r+1]]][t,] <- input[,r]
	}
#if (r == length(output)-1){
return(output)
#	}
}

# test code
# input <- matrix(seq(1,12),ncol=2,nrow=6)
#output <- list(id = 1, mat = matrix(0, ncol = 6, nrow = 2),
# 			mat2 = matrix(0, ncol = 6, nrow = 2))
# r <- 1:ncol(input)
# t <- 2
# list.write.fun(r,output,input,t)
# sapply(r,list.write.fun,output,input,1)



# +++++++++++++++++++++++++++++++++++++++++
# ------------------------------------
# Start of the run
#
# we are running over i timesteps and over j landscape cells
# GWthreshold = 0.10 m (set in gwt_input_parameters.txt)

big.fun <- function(N, stype, vtype, Rain, ETp, stream, gwheads, Zmean,
		rdir_in = "X:/vervoort/research/rcode/ecohydrology/2dmodelling",
		today.m = today) {
#print(gwheads)

# -----------------------------
# Preliminaries
# N: number of days for run
# stype: soiltype from soilfunction.r
# vtype: vector of vegetation types for vegfunction.r length no of cells
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
# --------------------------

# run read.fun
read.fun(rdir = rdir_in,stream.m = stream, gwheads.m = gwheads) # has no real output

  soilpar_in <- Soil(stype)
  n <- soilpar_in$n
  spec_y <- soilpar_in$spec_y


	# create a name for any output files
	NameRun<-paste(today.m,N,stype,sep="_")


  R <- Rain[,2] #vector with rain data
  deltat <- 1/12  
  DELT <- 1  
  finalN <- N #this will be overwritten when formula stops before N is reached  

    # introduce slopevector
	if(NY==1){  slopevector.m <- slope
    			}else{
    		slopevector.m<-matrixtovector(slope)  }

# -------- end preliminaries ----------------------
# xxxxxxxxxxxxxxxxxcheck
#  t<-1:(N*1/deltat) # needed??
#  ########
# D<-init_heads-bottom
# XXXXXXXXXXXXX check xxxxxxxxxxxxxxxxxxxx




# -------Create storage frames
	Storage_day <- list(id = 1:N, s = matrix(0,nrow=N,ncol=NX*NY), 
		qcap = matrix(0,nrow=N,ncol=NX*NY), Tgw = matrix(0,nrow=N,ncol=NX*NY),
		Tsoil = matrix(0,nrow=N,ncol=NX*NY), Ttotal = matrix(0,nrow=N,ncol=NX*NY), 
		Leakage = matrix(0,nrow=N,ncol=NX*NY), GWrech = matrix(0,nrow=N,ncol=NX*NY),
		static_stress = matrix(0,nrow=N,ncol=NX*NY), smloss = matrix(0,nrow=N,ncol=NX*NY),
		surfoff = matrix(0,nrow=N,ncol=NX*NY))
	Storage_Gridday <- data.frame(s = numeric(length=NX*NY), 
		qcap = numeric(length=NX*NY), Tgw = numeric(length=NX*NY),
		Tsoil = numeric(length=NX*NY), Ttotal = numeric(length=NX*NY), 
		Leakage = numeric(length=NX*NY), GWrech = numeric(length=NX*NY),
		static_stress = numeric(length=NX*NY), smloss = numeric(length=NX*NY),
		surfoff = numeric(length=NX*NY))

	Storage_Subday <- data.frame(s = numeric(length=1/deltat), 
		qcap = numeric(length=1/deltat), Tgw = numeric(length=1/deltat),
		Tsoil = numeric(length=1/deltat), Ttotal = numeric(length=1/deltat), 
		Leakage = numeric(length=1/deltat), GWrech = numeric(length=1/deltat),
		smloss = numeric(length=1/deltat),static_stress = numeric(length=1/deltat), 
		surfoff = numeric(length=1/deltat))


	GW_store <- list(id = 1:(N+1), heads = matrix(0,nrow=N+1,ncol=NX*NY), 
				GWdepths = matrix(0,nrow=N+1,ncol=NX*NY),
			 	delta_h = matrix(0,nrow=N+1,ncol=NX*NY),	
				qrivlat = matrix(0,nrow=N+1,ncol=NX*NY), 
				GWcum = matrix(0,nrow=N+1,ncol=NX*NY))
# ----------------------------------------------

# START OF THE LOOPS --------------------------#
# run a loop through the days
	for (t in 1:N) {
		if (t == 1) {
			# run gwt.exe function
			GW.out <- Gwt.fun(t, last_t_heads = gwheads, 
				RB = gwheads[NX], first=T, soilpar.m=soilpar_in, 
				GWcumvector = rep(0,NX), 
				slopevector = slopevector.m,
				DELT.m = DELT,
				NRBL.m = NRBL)  
			DELT <- GW.out$DELT
			# define initial soil moisture
			s_init <- rep(0.5,NX*NY)
		#	print(s_init) 
		} else {
			# check cumulative recharge function
			if (any(GW_store$GWcum[t-1,] >= GWthreshold)==T) {
			# if needed run gwt.exe function
				GW.out <- Gwt.fun(t,last_t_heads = GW_store$heads[t-1,], 
				first=F, RB = gwheads[NX], soilpar.m=soilpar_in, 
				GWcumvector=GW_store$GWcum[t-1,], 
				slopevector = slopevector.m,
				DELT.m = DELT,
				NRBL.m = NRBL)  
				# reset DELT
				DELT <- GW.out$DELT
				GW_store$GWcum[t,] <- 0
				print("gwt.exe is run")
			} # close recharge check loop
			# change DELT
			DELT <- GW.out$DELT <- DELT + 1
		} # close gw run loop
		# write away gw output 
		# if Gwt.fun wasn't run, it just writes the last one again
		r <- 1:(ncol(GW.out$out))
		#print(GW.out$out)
		GW_store <- list.write.fun(r, input=GW.out$out, 
			output = GW_store, t = t) 
		#str(GW_store)
		
	# Run a loop over the gridcells
    		for(j in 1:(NX*NY)) {
			# generate vegpar			
			vegpar <- Veg(vtype=vtype[j], soilpar=soilpar_in)
		  	# This is a variable based on Vervoort and van der Zee (2012)
     		vegpar$fs <- 0.5
        		vegpar$Ep <- ETp[t,2]
			# Define which model will be used (DR or no DR)
	  		if (vegpar$DR == TRUE) model_in <- "FB_new" else model_in <- rho_new_1 
			
			# run a loop for within day (this can be done for all x and y)
			for (p in 1:12) {		

				if (p == 1) { 
					# call rainfall function
					In <- Rain_in(p,R = R, n = soilpar_in$n, Zr = vegpar$Zr)
				} else {
					# rainfall = 0
					In <- 0
				}
		print(paste("p =",p))
	
				# run the water balance function
				s.old <-ifelse(p == 1,s_init[j],Storage_Subday$s[p-1])
				WB.out <- WB_fun(vegpar, In, fullday = t, 
					soilpar.m = soilpar_in,
					last_t_soil.sat = s.old, 
					Zmean.m = Zmean[j],
					deltat.m = deltat,
					model = model_in,
		   			Z.prev = ifelse(t==1,gwheads[j]*-100,GW_store$GWdepths[t-1,j]*-100),
		   			Z.m = ifelse(t==1,gwheads[j]*-100,GW_store$GWdepths[t,j]*-100))
				#Write water balance output to subdaily vector	
				Storage_Subday[p,] <- t(do.call(c,WB.out))
				# store s_init for gridcell till next day
				if (p == 12) s_init[j] <- Storage_Subday$s[p]
				} # close p loop
			# recalculate to daily output and write away
			Storage_Gridday[j,] <- apply(Storage_Subday,2,sum)	
			Storage_Gridday[j,c(1,9)] <- apply(Storage_Subday,2,mean)[c(1,9)]	
			# reset Storage_subday??
print(paste("j =",j))
		} # close j loop
#		print(s_init)
 		# daily values across grid
		r <- 1:ncol(Storage_Gridday)
		Storage_day <- list.write.fun(r, input=Storage_Gridday, 
		 		output = Storage_day, t=t)
		#str(Storage_day)
		# Accumulate the GW recharge, this is now in GW_store
		GW_store$GWcum[t,] <- GW_store$GWcum[t,] + Storage_day$GWrech[t,]*0.01 # m/day		
		print(paste("t =",t))
	} # close t loop

# ------------ END OF LOOPS ---------------------

# Run gwt_fun one last time
	GW.out <- Gwt.fun(t = N+1, last_t_heads = GW_store$heads[t-1,], first=F,
				RB = gwheads[NX], soilpar.m=soilpar_in, 
				GWcumvector=GW_store$GWcum[t-1,], 
				slopevector = slopevector.m,
				DELT.m = DELT,
				NRBL.m = NRBL)
	# write away gw output 
	r <- 1:(ncol(GW.out$out))
	GW_store <- list.write.fun(r, input=GW.out$out, 
	 		output = GW_store, t = N + 1) 

# ---------------------------------------------
# Write and store output
# -------------------------------------
# write away values to some files (WB etc)
Storage.out <- do.call(cbind,Storage_day)
write.csv(Storage.out,paste(NameRun,"DailyOutput_table.csv",sep="_"),
		row.names=FALSE)
GWstore.out <- do.call(cbind,GW_store)
write.csv(GWstore.out,paste(NameRun,"GWOutput_table.csv",sep="_"),
		row.names=FALSE)

# Create a water balance table
# Not that important
Wb <- apply(Storage.out,2,sum)
Wb[c(1,8)] <- apply(Storage.out,2,mean)[c(1,8)]
#print(Wb)

write.table(Wb,paste(NameRun,"Output_WB_table.csv",sep="_"),row.names=FALSE,
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
