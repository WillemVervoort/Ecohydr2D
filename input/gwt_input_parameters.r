# ----------------
# Grid
# ---------------
NX <- 42
NY <- 1
NRBL <- 2    #number of river blocks
DELX <- c(20,rep(10,10), rep(25,10), rep(50,10),rep(500,10),20)
DELY <- 100   
ITIM <- 1
distancetoriver=c(rep(0,NRBL),cumsum(DELX[NRBL:length(DELX)]))
# Define slope grid
dslope_x <- 0 # m difference in x-direction
dslope_y <- 0 # m difference in y-direction
# define the slopw
slope <- slopefun(NX,NY,dslope_x,dslope_y)
# ----------------------------

# --------------------
# Groundwater
# --------------------
# Depth of water table:
init_heads <- gwheads #meters
bottom <- -25   #meters
bottommatrix<-matrix(bottom,NX,NY)
bottomvector<-matrixtovector(bottommatrix)


GWthreshold<- 0.0 #m 
DELTcrit<- 21   #days
# -------------------


########################
#river: 
########################
criver <- c(RES,0.1) #resistance (m) # WV says Increase this to limit leakage from river
Ariver <- rep(100,NRBL) # river area (m^2) # this is a guess
riverheads <- (stream[,2] - 3)  # water head river (m). 
# This assumes no overbank flows in the period. 
#creates riverheads vector

if(NRBL==0){hriver <- 0} # should have some value, otherwise: formula breaks off

#alternative:river on all outside blocks. NOT IMPORTANT
#boundaryriver <- function(NX, NY, criver, Ariver) {
#  NRBL<- 2*NX+ 2*(NY-2)
#  rivermatrix <- matrix(0,NRBL,4) 
#  rivermatrix[,3] <- Ariver
#  rivermatrix[,4] <- criver
#  i<-1
#  for(j in 1:NY){
#    rivermatrix[i,1] <- 1
#    rivermatrix[i,2] <- j
#    i<-i+1
#    rivermatrix[i,1] <- NX
#    rivermatrix[i,2] <- j
#    i<-i+1
#  }
#  for (j in 2:(NX-1)){
#    rivermatrix[i,1] <- j
#    rivermatrix[i,2] <- 1
#    i<-i+1
#    rivermatrix[i,1] <- j
#    rivermatrix[i,2] <- NY
#    i<-i+1    
#  }    
#  return(rivermatrix)     
#}


#rivermatrix <- boundaryriver(NX=NX, NY=NY, criver=criver, Ariver=Ariver)

#########################
#RAIN (Older stuff from Joep, to model rainfall)
########################

# Different way of writing the Rainfall function
#Precip <- function(time,alpha,lambda,delta) {
#       # generate a vector of times between rainfall events > delta
#       f_P<-round(rexp(time,lambda*exp(-delta/alpha))) # vector of times between rainfall occurrences (equation 4 & 8)
#       # generate a binary vector from this (impulse function)
#       binary.vec <- unlist(lapply(1:time,function(i) c(rep(0,f_P[i]),1)))
#       R <- rexp(length(binary.vec),1/alpha)*binary.vec 
#       return(R[1:time])
#    }
#    # Delta=interceptie: functie van vegetatie. 
#    BourkeGrass <- Precip(time=15000, alpha=0.74, lambda=0.14, delta=0.1)
#    BourkeTrees <- Precip(time=15000, alpha=0.74, lambda=0.14, delta=0.2)
#    MoreeGrass <- Precip(time=15000, alpha=0.89, lambda=0.21, delta=0.1)
#    MoreeTrees <- Precip(time=15000, alpha=0.89, lambda=0.21, delta=0.2)
#    KatherineGrass <- Precip(time=15000, alpha=1.47, lambda=0.22, delta=0.1)
#    KatherineTrees <- Precip(time=15000, alpha=1.47, lambda=0.22, delta=0.2)
#    
#    #Rain[1:NRain,2]<-Precip(time=NRain, alpha=0.74, lambda=.14, delta=vegpar$delta)
#    
#    Norain<-function(N){
#      Rain[1:N,2]<-0
#      print(Rain[1:N,])
#      }
#    #Rain<-Norain(NRain)
#
# Monthly solar exposure at Goodooga, NSW in MJ/m2/day
#Rs <- 10^6*cbind(1:12,c(27.3,24.8,22.5,17.7,13.9,11.8,12.9,16.2,20.4,23.7,25.7,27.5))
#e_star <- 0.22

# Expand the Im vector to same length as rainfall
#temp <- cbind(format(as.Date(Rain_data[,1],"%d/%m/%Y"),"%m"),rep(0,nrow(Rain_data)))
#for (i in 1:12) {
#  temp[as.numeric(temp[,1])==i,2] <- (0.5/e_star)*Rs[i,2]
#}
#I_m <- as.numeric(temp[1:nrow(Rain),2])



