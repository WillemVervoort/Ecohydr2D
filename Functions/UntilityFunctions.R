# utility functions

# Do a stack operation, Joep actually wrote a function for this
stackfun <- function(data,NX) { # this needs to be adapted relative to what you want
  df <- data.frame(days=rep(data$id,NX), Rain=rep(data$Rain,NX),
                   loc=rep(data$x,each=length(data$id)),
                   River=rep(data$Stream,NX), 
                   Ts = stack(as.data.frame(data$Tsoil[,1:NX]))[,1],
                   Tg = stack(as.data.frame(data$Tgw[,1:NX]))[,1],
                   Ttotal = stack(as.data.frame(data$Ttotal[,1:NX]))[,1], 
                   s = stack(as.data.frame(data$s[,1:NX]))[,1],
                   qcap = stack(as.data.frame(data$qcap[,1:NX]))[,1],
                   gwlevel = stack(as.data.frame(data$gwlevel[,1:NX]))[,1],
                   gwrech = stack(as.data.frame(data$GWrech[,1:NX]))[,1],
                   stress = stack(as.data.frame(data$static_stress[,1:NX]))[,1])
  return(df)
}

# Read in function for sourcing other scripts
read.fun <- function(rdir = rdir_in) {
  source(paste(rdir,"soilfunction.r",sep="/"))
  source(paste(rdir,"vegfunction.r",sep="/"))
}
read.fun1 <- function(rdir = rdir_in, stream.m = stream, gwheads.m = gwheads,
                      Res = NULL) {
  attach(list(gwheads=gwheads.m,stream=stream.m, RES= Res))
  
  source(paste(rdir,"gwt_input_parameters.r",sep="/")) #You need to adjust values in here as well!!
  source(paste(rdir,"Fluxfunctions2Dmodel.R",sep="/"))  #
  source(paste(rdir,"define_input.r",sep="/"))      #also the parameters are called
  
  # detach again to make sure there is no confusion
  detach(list(gwheads = gwheads.m, stream = stream.m, RES = Res))
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
# t <- 1
# list.write.fun(r,output,input,t)
# sapply(r,list.write.fun,output,input,1)
