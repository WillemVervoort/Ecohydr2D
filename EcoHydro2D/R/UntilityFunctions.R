# utility functions

# Do a stack operation, Joep van der Zanden actually wrote a function for this
stackfun <- function(data,NX) { 
   # this needs to be adapted relative to what you want
    #browser()
    df <- data.frame(days=rep(data$id,NX), rain=rep(data$rain,NX),
                   location=rep(data$x[(length(data$x)-NX+1):length(data$x)],each=length(data$id)),
                   vtype=rep(data$vtype[(length(data$vtype)-NX+1):length(data$vtype)],each=length(data$id)),
                  stream=rep(data$stream,NX),qrivlat=rep(data$qrivlat,NX),  
                   Ts = stack(as.data.frame(data$Tsoil[,1:NX]))[,1],
                   Tg = stack(as.data.frame(data$Tgw[,1:NX]))[,1],
                   Ttotal = stack(as.data.frame(data$Ttotal[,1:NX]))[,1], 
                   s = stack(as.data.frame(data$s[,1:NX]))[,1],
                   qcap = stack(as.data.frame(data$qcap[,1:NX]))[,1],
                   gwlevel = stack(as.data.frame(data$gwlevel[,1:NX]))[,1],
                   gwrech = stack(as.data.frame(data$gwrech[,1:NX]))[,1],
                   leakage = stack(as.data.frame(data$leakage[,1:NX]))[,1],
                   surfoff = stack(as.data.frame(data$static_stress[,1:NX]))[,1],
                   smloss = stack(as.data.frame(data$static_stress[,1:NX]))[,1],
                   stress = stack(as.data.frame(data$static_stress[,1:NX]))[,1])
return(df)
}

# # Read in function for sourcing other scripts
# this needs to be rewritten to be based on ecohydro2d.options
read.fun1 <- function(stream.m, gwheads.m) {
  NX <- ecohydro2d.options()$NX
  NY <- ecohydro2d.options()$NY
  DELX <- ecohydro2d.options()$DELX
  NRBL <- ecohydro2d.options()$NRBL
  dslope_x <- ecohydro2d.options()$dslope_x
  dslope_y <- ecohydro2d.options()$dslope_y
  RES <- ecohydro2d.options()$RES
  #browser()
  bottommatrix<-matrix(ecohydro2d.options()$bottom,NX,NY)
  bottomvector<-matrixtovector(bottommatrix)
  # This assumes no overbank flows in the period. 
  # creates riverheads vector
  if(NRBL==0) hriver <- 0 else hriver = NULL # should have some value, otherwise: formula breaks off
  distancetoriver=c(0,cumsum(DELX[1:length(DELX)]))
  slope <- slopefun(NX,NY,
            dslope_x,dslope_y)
  
  # Depth of water table:
  riverheads <- stream.m[,2]  # water head river (m).
  init_heads <- gwheads.m # groundwater heads (m).
  # river parameters
  criver = c(RES,0.1) #resistance (m) # WV says Increase this to limit leakage from river
  Ariver <- rep(100,NRBL) # river area (m^2) # this is a guess
  
  return(list(bottommatrix=bottommatrix,bottomvector=bottomvector,
                       slope=slope, riverheads=riverheads,
                       distancetoriver = distancetoriver,
                       init_heads = init_heads, criver = criver,
                       Ariver = Ariver, hriver=hriver))
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

list.write.fun2 <- function(output,input,t) {
  
  # input and output are a dataframe and a list
  # i is the time counter
  # r is a counter along the columns of inputs
  for (r in 1:nrow(input)) {
    output[[names(output)[r+1]]][t,] <- input[r,]
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
# t <- 1
# list.write.fun(r,output,input,t)
# sapply(r,list.write.fun,output,input,1)

## Write important data to file
SaveParameterData <- function(name){
  Out<- cbind(NX,NY,NRBL,t(DELX),DELY,init_heads,bottom,dslope_x, dslope_y,
              criver, Ariver, riverheads)
  write.table(Out,name,row.names=FALSE,col.names=TRUE,sep=",")
}

matrixtovector<-function(inputmatrix){
  vect<-vector(mode="numeric", length=(nrow(inputmatrix)*ncol(inputmatrix)))
  a<-1
  for (i in 1:nrow(inputmatrix)){
    for (j in 1:ncol(inputmatrix)){
      vect[a]<-inputmatrix[i,j]
      a<-a+1
    }
  }
  return(vect)
}


vectortomatrix<-function(inputvector, Ncol){
  Nrow<-length(inputvector)/Ncol
  outputmatrix <- matrix(0,nrow=Nrow,ncol=Ncol)
  a<-1
  for (j in 1:ncol(outputmatrix)){
    for (i in 1:nrow(outputmatrix)){
      outputmatrix[i,j] <-  inputvector[a]
      a<-a+1
    }
  }
  return(outputmatrix)
}



slopefun <- function(NX,NY,dslope_x,dslope_y) {
  slope<-matrix(0,NY,NX)
  if(NY==1){slope<-seq(0,dslope_x,dslope_x/(NX-1))
  } else {
    for(i in 1:NX){
      for(j in 1:NY){
        slope[j,i]=(i-1)*(dslope_x)/(NX-1)+(j-1)*(dslope_y)/(NY-1)
      }
    }
  }
  return(slope)
}


