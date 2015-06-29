today <- format(Sys.Date(),"%Y%m%d")
setwd("E:/Back-up Mijn Documenten/Mijn documenten 2010-09-10 10;07;24/Universiteit/Thesis/My Documents/Thesis/Working directory")
#rdir <- "E:/Back-up Mijn Documenten/Mijn documenten 2010-09-10 10;07;24/Universiteit/Thesis/My Documents/Thesis/Working directory"
rdir <- "D:/Joep/Thesis/Working directory"
#source(paste(rdir,"gwt_input_parameters.txt",sep="/")) # 

# In advance, define all parameters in gwt_input_parameters.txt

      
#Model: function of (N)
# ALL ONE BIG FUNCTION DOWN TO LINE 413
Coupling <- function(N,stype,vtype,climatetype, MeanRain) {
#read initial values and call other codes

source(paste(rdir,"soilfunction.r",sep="/"))
source(paste(rdir,"vegfunction.r",sep="/"))


soilpar <- Soil(stype)
vegpar <- Veg(vtype=vtype, soilpar=soilpar)
DR <- vegpar$DR

# SETS THE DEPTH OF WATER TABLE (CONTROLLED BY DEEP)
if(deep==FALSE){riverheads <-  vegpar$Zr/-100 - 2}else{
  riverheads <-  vegpar$Zr/-100 - 3.5}  

source(paste(rdir,"gwt_input_parameters.txt",sep="/")) #out because otherwise: overwrite veg,soil,climate data 
#source(paste(rdir,"climatefunction.r",sep="/"))
source(paste(rdir,"flux_equations_NEW.r",sep="/"))  #!!!!! NEW means: adjusted formula.
source(paste(rdir,"define_input.r",sep="/"))      #also the parameters are called

NameRun<-paste(today,N,stype,vtype,climatetype, MeanRain)
#!!!!!!!!!!!!!!!! Temporary adjusting the initial values!!
#if(stype=="L Med Clay" && vtype=="Grass" && climatetype =="Moree"){ init_heads <- -9.10}
#if(stype=="L Med Clay" && vtype=="TreesDR" && climatetype =="Moree"){ init_heads <- -59.00}
#if(stype=="L Med Clay" && vtype=="TreesNoDR" && climatetype =="Moree"){ init_heads <- -9.80}
#if(stype=="Loamy Sand" && vtype=="Grass" && climatetype =="Moree"){ init_heads <- -10}
#if(stype=="Loamy Sand" && vtype=="TreesDR" && climatetype =="Moree"){ init_heads <- -86}
#if(stype=="Loamy Sand" && vtype=="TreesNoDR" && climatetype =="Moree"){ init_heads <- -10.4}
#if(stype=="S Clay Loam" && vtype=="Grass" && climatetype =="Moree"){ init_heads <- -14.3}
#if(stype=="S Clay Loam" && vtype=="TreesDR" && climatetype =="Moree"){ init_heads <- -63}
#if(stype=="S Clay Loam" && vtype=="TreesNoDR" && climatetype =="Moree"){ init_heads <- -14.9}
#if(stype=="L Med Clay" && vtype=="Grass" && climatetype =="Katherine"){ init_heads <- -2.60}
#if(stype=="Loamy Sand" && vtype=="Grass" && climatetype =="Katherine"){ init_heads <- -2.8}
#if(stype=="S Clay Loam" && vtype=="Grass" && climatetype =="Katherine"){ init_heads <- -3.7}
#if(stype=="Loamy Sand" && vtype=="TreesDR" && climatetype =="Katherine"){ init_heads <- -2.50}

#### Rain
Rain_data<-read.csv(paste(rdir,"Rainfall data",paste(climatetype,vtype,".csv",sep=""),sep="/"))
Rain <- data.frame(dates=as.Date(Rain_data[,1], "%d/%m/%Y"), Rain=Rain_data[,2]) # already in cm/d
if(MeanRain==TRUE) Rain[,2]<-mean(Rain[,2])
temp <- cbind(format(as.Date(Rain_data[,1],"%d/%m/%Y"),"%m"),rep(0,nrow(Rain_data)))
# add reading in the potential ET data

# for (i in 1:12) {
#   temp[as.numeric(temp[,1])==i,2] <- (0.5/e_star)*Rs[i,2]
# }
# I_m <- as.numeric(temp[1:nrow(Rain),2])


riverheads <-  vegpar$Zr/-100 - 2
heads<-matrix(0,nrow=(N+1),ncol=(NX*NY)) 

# initial heads
heads[1,]<-riverheads
 
if(length(DELX)>1){
  est_leakage<-0.0*sum(Rain[,2])/25000  # start with analytical approximation. If 0, a flat gw table is used
  total_xx <- 2*sum(DELX[2:NX])
  for(ii in 2:length(DELX)){                 #est_leakage in cm/d. K_s in cm/d
    heads[1,ii]<-sqrt(est_leakage/soilpar$K_s*(total_xx*distancetoriver[ii]-distancetoriver[ii]^2) +
      (heads[1,ii]-bottom)^2)+bottom
    }
  sum(DELX)
}

  #first define initial values of BigFooFunction. Original vectors become matrices (1D to 2D)
  if (DR == TRUE) model <- "Fun_FB" else model <- rho_new_1 # "Fun_E" #
  sinit<-vegpar$s_star
  Zr <- vegpar$Zr
  n <- soilpar$n
  spec_y <- soilpar$spec_y
  R <- Rain[(1:N),2] #vector with rain data
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
  LAI <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  DM <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  Aact <- matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  static_stress<-matrix(0,nrow=N*1/deltat,ncol=(NX*NY))
  ss[1,] <- sinit  # You need a begin value
  static_stress[1,]<- 0
  Tg[1,] <- 0
  Ts[1,] <- 0
  Leakage[1,]<-0
  smloss[1,] <- 0 
  LAI[1,] <- vegpar$LAI
  DM[1,] <- 0.5
  Tt[1,] <- (1-exp(-0.45*LAI[1,]))*vegpar$Ep
  qcap[1,] <- 0
  DORM <- FALSE
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
  GWdepths[1,]<-heads[1,]-t(slopevector)

#  # create table with initial values. 
#  initial_table<- data.frame(days = 0, heads=t(as.numeric(heads[1,])), GWdepth=t(as.numeric(GWdepths[1,])),
#      Rain=0, smloss=t(as.numeric(smloss[1,])), Ts=t(as.numeric(Ts[1,])),Tg=t(as.numeric(Tg[1,])), 
#      qcap=t(as.numeric(qcap[1,])), GWrech=t(as.numeric(GWrech[1,])),
#      surfoff=t(as.numeric(surfoff[1,])),s=t(as.numeric(ss[1,])))
#          # LAI=t(as.numeric(LAI[1,])),DM=t(as.numeric(DM[1,])),  Aact=t(as.numeric(Aact[1,])), 

#Define parameters in Toons model{
  #adjust parameters here if necessary: ksat, grid, recharge, ...
  #create grid. start with zero RECHarge vector for first run Toons model.
  define_gwt_gridinput(NX, NY, NRBL, widthxvector, widthyvector)#, NB, XB1, YB1, XB2, YB2, XB3, YB3, XB4, YB4)
  define_gwt_timestepinput(ITIM, DELT, RECH=GWrech[1,], hriver=riverheads)
  define_gwt_hydroinput(heads[1,], bottomvector,ksat=soilpar$K_s*0.01,spec_y=soilpar$spec_y)  #K_s from cm/d to m/d
  if(NRBL>0) {define_gwt_riverinput(rivermatrix=rivermatrix)}
  #Run Toons model once to get the proper head output / the first Z input
#  system(paste("'",paste(rdir,"GWT.exe",sep="/"),"'"))
  system("GWT.exe")

  #Loop: time
  i<-1
  STOP<-FALSE  
  while(i<=(nrow(ss)-1)&STOP==FALSE) {  
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
     define_gwt_timestepinput(ITIM=(i*deltat+1), DELT=DELT, RECH=t(GWcumvector)/DELT, hriver=riverheads)
     system("GWT.exe")
     DELT<-1
     heads[floor(i*deltat)+1,]<-matrixtovector(read.table(paste(rdir,"gwt_headoutput",sep="/"),fill=TRUE))[1:(NX*NY)] #Fill: read also last lines when table is no square matrix
     GWdepths[floor(i*deltat)+1,]<-heads[floor(i*deltat)+1,]-t(slopevector)
     delta_h[floor(i*deltat)+1,]<- heads[floor(i*deltat)+1,]-heads[floor(i*deltat),]
     qrivlat[floor(i*deltat)+1,]<-(delta_h[floor(i*deltat)+1,]-GWcumvector/spec_y)*100*spec_y
     GWcummatrix[floor(i*deltat)+1,] <- t(GWcumvector)*100
     GWcumvector[1:(NX*NY)]<-0 #Reset GWcumvector
     } 
     if (((floor(i*deltat)+1)!=(i*deltat+1))) {  #Means: every time step except the 12th  
     In <- 0
     }
     ## Include stop code: 
      if(max(heads[floor(i*deltat)+1,])>-Zr/100) {
      finalN <- i*deltat
      print(paste("GW reaches root zone for time step = ",finalN))
      print(c(stype,climatetype,vtype))
      
#      heads[N+1,]<-(-Zr/100)
      STOP<-TRUE
      }
  # calculate root zone model for each Z in Toons model --> Loop: space
     if(STOP==FALSE){
     for(j in 1:(NX*NY)) {
        # define Z values 
        Z <- GWdepths[floor(i*deltat)+1,j]*-100
        # The Big Foo Loop pasted!
        #  Empty[i]<-n*Zr*(1-ss[i])

        Phi[i,j]<-In
 #       print(Phi[i,j])
 #       print(ss[i,j])
        if(Phi[i,j]>(1-ss[i,j])){
          Phi[i,j]<- 1-ss[i,j]
          surfoff[i+1,j]<- (In-Phi[i,j])*n*Zr
          }
      # in here read in the actual pot ET into vegpar$Ep
        # vegapr$Ep <- ETdata[i]
        loss <- do.call(model,list(s=ss[i,j],soilpar=soilpar,vegpar=vegpar,Z=Z,Z.prev=(GWdepths[max(1,floor(i*deltat)),j]*-100))) #i should be at least 1
#        ss[i+1,j]<-ss[i,j]+Phi[i,j]-loss[1]/(n*Zr)*deltat
        
        # update the water balance

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
#        qcap[i+1,j] <- loss[4]#/(n*Zr)*deltat           ###!!!
        Ts[i+1,j] <- loss[2]#/(n*Zr)*deltat
        Tt[i+1,j] <- Tg[i+1,j] + Ts[i+1,j]
        Leakage[i+1,j] <- ifelse(loss[5]>0,loss[5],0)

        smloss[i+1,j] <- loss[1]  - (qcap[i+1,j] - loss[4]) #increase with difference between originally calculated and new (pot.reduced) capillary flux
        GWrech[i+1,j] <- Leakage[i+1,j]-loss[3]-qcap[i+1,j]

        # Now grow some biomass and LAI
        # Calculate FPAR (Lu et al. 2001)
        # Assume extinction coefficient of 0.4 (SWAP)
#        FPAR <- vegpar$fc.coef*100*(1-exp(-0.45*LAI[i,j]))
#        #print(paste("FPAR =",FPAR))
#        # Calculate Potential assimilation (Haxeltine et al. 1996)
#        sigm <- (1-0.08*deltat)^0.5
#        Ap <- FPAR*sigm*vegpar$veg.psi*I_m[floor(i*deltat)+1]*deltat
#        # And Actual Assimilation using ET
#        Aact[i+1,j] <- Tt[i+1,j]/((1-exp(-0.45*LAI[i,j]))*vegpar$Ep)*Ap*deltat
#        # Subtract respiration
#        A_act <- (Aact[i+1,j]-0.4*Aact[i+1,j])
#        # Convert to dry Matter [g] :
#        DM[i+1,j] <- 2*A_act/10^6
#        # Dry Matter in the foliage [g / m^2 ground] :
#        stress <- ifelse(ss[i,j] < vegpar$s_star,(vegpar$s_star-ss[i,j])/(vegpar$s_star-vegpar$s_w),0)
#        if (vegpar$dorm==TRUE && vegpar$LAI < vegpar$dorm.cutoff*LAI[1,j]) {
#          DMf <- DM[i+1,j]* vegpar$AboveGround_PF * vegpar$Leaf_PF
#          # Estimate of LAI from the Dry Matter content in the foliage :
#          # SLA = Specific Leaf Area [m^2 leaf / kg DM]
#          #Dead.leaf[i] <- vegpar$Mort*deltat*LAI[i]/vegpar$SLA
#          dLAI <- vegpar$SLA * DMf/1000  - vegpar$Mort*deltat*LAI[i,j]
#          DORM <- TRUE
#        } else {
#          DMf <- DM[i+1,j]* vegpar$AboveGround_PF * vegpar$Leaf_PF
#          # Estimate of LAI from the Dry Matter content in the foliage :
#          # SLA = Specific Leaf Area [m^2 leaf / kg DM]
#          dLAI <- vegpar$SLA * DMf/1000  - vegpar$Mort*(1+stress)*deltat*LAI[i,j]
#        }
#        if (stress == 0 && DORM == TRUE) {
#          DMf <- DM[i+1,j]* vegpar$AboveGround_PF # lots of growth
#          dLAI <- vegpar$SLA * DMf/1000 #- vegpar$Mort*deltat*LAI[i]
#        }   
#        LAI[i+1,j] <- LAI[i,j] + dLAI
#        if (LAI[i+1,j] == LAI[1,j] && DORM == TRUE) DORM <- FALSE
#        vegpar$LAI <- LAI[i+1,j]
     GWcumvector[j] <-GWcumvector[j] + GWrech[i+1,j]*deltat*0.01 #so in m/d
     }
  
 #     print(heads[i,])
  i<-i+1
 }}
# execute GWT model one last time. + create final heads
   define_gwt_timestepinput(ITIM=(N+1), DELT=DELT, RECH=t(GWcumvector)/DELT, hriver=riverheads)
   system("GWT.exe")     
   heads[N+1,]<-matrixtovector(read.table(paste(rdir,"gwt_headoutput",sep="/"),fill=TRUE))[1:(NX*NY)] #Fill: read also last lines when table is no square matrix
   GWdepths[N+1,]<-heads[N+1,]-t(slopevector) #heads & GWdepths N+1 are not in output table. Just used to calculate dh for waterbalance.
   GWcummatrix[N+1,]<-t(GWcumvector)*100
   delta_h[N+1,]<- heads[N+1,]-heads[N,]    
   qrivlat[N+1,]<-(delta_h[N+1,]-GWcumvector/spec_y)*100*spec_y

                       

  # create output data
        days <- sort(rep(1:N,1/deltat))
        sat_out <- aggregate(ss,list(day=days),mean)
        static_stress_out<-aggregate(static_stress,list(day=days),mean)
        Ts_out <- aggregate(Ts,list(day=days),mean)
        Tg_out <- aggregate(Tg,list(day=days),mean)
        smloss_out <- aggregate(smloss,list(day=days),mean)
        surfoff_out <- aggregate(surfoff,list(day=days),mean)
        qcap_out <- aggregate(qcap,list(day=days),mean)
        leak_out <- aggregate(Leakage,list(day=days),mean)
        GWrech_out <- aggregate(GWrech,list(day=days),mean)  # back to cm/d
        LAI.out <- aggregate(LAI,list(day=days),mean)
        DM.out <- aggregate(DM,list(day=days),mean)
        Aact.out <- aggregate(Aact,list(day=days),mean)         

# write important values to hard disk (s, Z, fluxes Q, ET)
#create output_table, start with line with initial values
  #write.table(initial_table,paste(rdir,"Output_table",sep="/"),row.names=FALSE,col.names=TRUE,sep=" ")
#create table with output results       
output_table<- data.frame(days=sat_out[,1], Heads=heads[1:N,1:(NX*NY)], GWdepth=GWdepths[1:N,1:(NX*NY)],
    Rain=R, smloss=smloss_out[,2:(NX*NY+1)],Ts=Ts_out[,2:(NX*NY+1)],  Tg=Tg_out[,2:(NX*NY+1)], 
    qcap=qcap_out[,2:(NX*NY+1)], leak=leak_out[,2:(NX*NY+1)], GWrech=GWrech_out[,2:(NX*NY+1)],
    GWrechcum=GWcummatrix[2:(N+1),1:(NX*NY)], deltaH=delta_h[2:(N+1),1:(NX*NY)],surfoff=surfoff_out[,2:(NX*NY+1)], 
    qrivlat=qrivlat[2:(N+1),1:(NX*NY)], s=sat_out[,2:(NX*NY+1)], st_stress=static_stress_out[,2:(NX*NY+1)])
#  LAI=LAI.out[,2:(NX*NY+1)], DM=DM.out[,2:(NX*NY+1)], Aact=Aact.out[,2:(NX*NY+1)],
#paste output results in output_table
write.csv(output_table,paste(rdir,paste("Output_table",NameRun,sep=""),sep="/"),row.names=FALSE)
     

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
## get rid of the first 2000 days for long runs (N>20000)
if(N>20000) {
sum_rain<-sum(R[2001:finalN])
for(i in 1:(NX*NY)){
  sum_smloss[i]<-sum(smloss_out[2001:finalN,i+1])
  sum_Ts[i]<-sum(Ts_out[2001:finalN,i+1])
  sum_Tg[i]<-sum(Tg_out[2001:finalN,i+1])
  sum_qcap[i]<-sum(qcap_out[2001:finalN,i+1])
  sum_leak[i]<-sum(leak_out[2001:finalN,i+1])
  sum_GWrech[i]<-sum(GWrech_out[2001:finalN,i+1]) 
  sum_surfoff[i]<-sum(surfoff_out[2001:finalN,i+1])
  sum_delta_h[i]<-sum(delta_h[1:(N+1),i])
  deltas_real[i]<- sat_out[finalN,i+1]-ss[2001,i]
  deltas_calc[i]<-(sum_rain-sum_smloss[i]-sum_surfoff[i])/(vegpar$Zr*n)
  mean_s.T[i] <- mean(sat_out[2001:finalN,i+1])
  mean_s.2[i] <- mean(sat_out[(2000+(floor((finalN-2000)/2))):finalN,i+1])
  mean_Z.T[i] <- mean(heads[2001:(finalN+1),i])            # calculate mean_Z
  mean_Z.2[i] <- mean(heads[(2000+(floor((finalN-2000)/2))):(N+1),i]) # calculate mean_Z of second half of run - check if the mean changes.
  var_s.T[i] <- var(sat_out[2001:finalN,i+1])            # calculate var_s
  var_s.2[i] <- var(sat_out[(2000+(floor((finalN-2000)/2))):finalN,i+1]) # calculate var_s of second half of run - check if the mean changes.
  var_Z.T[i] <- var(heads[2001:(N+1),i])            # calculate var_Z
  var_Z.2[i] <- var(heads[(2000+(floor((finalN-2000)/2))):(N+1),i]) # calculate var_Z of second half of run - check if the mean changes.
  cell[i]<-i
  sum_qrivlat[i]<-sum(qrivlat[2001:(N+1),i])             
 }
 }


WB_table<- data.frame(cell=cell, slope=slopevector, deltah=sum_delta_h, 
    sum_rain=sum_rain, sum_smloss=sum_smloss, sum_Ts=sum_Ts, sum_Tg=sum_Tg, sum_qcap=sum_qcap, 
    sum_leak=sum_leak, sum_GWrech=sum_GWrech, sum_qrivlat=sum_qrivlat, sum_surfoff=sum_surfoff, 
    deltas_real=deltas_real, deltas_calc=deltas_calc, mean_s.T=mean_s.T, mean_s.2=mean_s.2,
    var_s.T=var_s.T, var_s.2=var_s.2, mean_Z.T=mean_Z.T, mean_Z.2=mean_Z.2, var_Z.T=var_Z.T,
    var_Z.2=var_Z.2) 
write.table(WB_table,paste(rdir,paste("Output_WB_table",NameRun,sep=""),sep="/"),row.names=FALSE,
    col.names=TRUE,sep=" ")

## Output is values for the first cell!
#print(data.frame(days=sat_out[,1],Heads=heads[1:N,1],GWdepth=GWdepths[1:N,1], Rain=R, 
#    smloss=as.numeric(smloss_out[,2]),Ts=as.numeric(Ts_out[,2]),Tg=as.numeric(Tg_out[,2]), 
#    qcap=as.numeric(qcap_out[,2]),GWrech=as.numeric(GWrech_out[,2]),GWrechcum=as.numeric(GWcummatrix[2:(N+1),1]),
#    RivLat=as.numeric(qrivlat[2:(N+1),1]), deltaH=delta_h[2:(N+1),1],surfoff=surfoff_out[,2],
#    s=as.numeric(sat_out[,2]), st_stress=as.numeric(static_stress_out[,2]))) 
##        LAI=as.numeric(LAI.out[,2]),DM=as.numeric(DM.out[,2]), Aact=as.numeric(Aact.out[,2]),    
##print(heads[N+1,])

#Output for stack formula:
return(data.frame(x=distancetoriver,heads0=heads[1,],headshw=heads[floor(N/2),],
    headsN=heads[N+1,],meanZ=mean_Z.T,varZ=var_Z.T,mean_s=mean_s.T,var_s=var_s.T,
    meanZ2=mean_Z.2,varZ2=var_Z.2,mean_s2=mean_s.2,var_s2=var_s.2,
    VEG=vtype,SOIL=stype,K_s=soilpar$K_s,CLIMATE=climatetype, MeanRain=MeanRain))
}

# Now run a few different simulations
# THIS SETS THE DEPTH OF THE INITIAL WATER TABLE
deep <- FALSE

# TO MODEL DEEP ROOTS WITH TREES YOU NEED TO CHANGE vtype to
# vtype = "TreesDR" otherwise "TreesNoDR"
Coupling(N=10,stype="L Med Clay", vtype="Grass",climatetype="clim1",MeanRain=FALSE)

# assign to look at output
result <- Coupling(N=10,stype="L Med Clay", vtype="Grass",climatetype="clim1",MeanRain=FALSE)
result$x # should give you distance to river


########## make a plot
#par(mfrow=c(2,1),mar=c(4,4,2,5)) # Leave space for z axis
#plot(Out$days,Out$Rain, xlab = "Time (days)", ylab="Rain (cm/d)")
#par(new = TRUE)
#plot(Out$days,Out$Heads, type = "l", axes=FALSE, bty = "n", xlab = "", ylab="")
#axis(side=4, at = pretty(range(Out$Heads)))
#mtext("Groundwater depth (m)", side=4, line=3)
#par(mar=c(4,4,2,5))
#plot(Out$days,Out$s, type="l", col="blue", xlab = "Time (days)", ylab="s (-)")
#par(new = TRUE)
#plot(Out$days,Out$Heads, type = "l", axes=FALSE, bty = "n", xlab = "", ylab="")
#axis(side=4, at = pretty(range(Out$Heads)))
#mtext("Groundwater depth (m)", side=4, line=3)
