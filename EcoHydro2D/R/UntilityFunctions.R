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

# # Read in function for sourcing other scripts
# read.fun <- function(rdir = rdir_in) {
#   source(paste(rdir,"soilfunction.r",sep="/"))
#   source(paste(rdir,"vegfunction.r",sep="/"))
# }
read.fun1 <- function(input_dir, stream.m = stream, gwheads.m = gwheads,
                      Res = NULL) {
  attach(list(gwheads=gwheads.m,stream=stream.m, RES= Res))
  
  source(paste(input_dir,"gwt_input_parameters.r",sep="/")) #You need to adjust values in here as well!!
#  source(paste(rdir,"define_input.r",sep="/"))      #also the parameters are called
  
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
  Out<- cbind(NX,NY,NRBL,t(DELX),DELY,init_heads,bottom,GWthreshold,DELTcrit,dslope_x, dslope_y,
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

# 
# if(NX==1){
#   widthxvector<-DELX
#   distancetoriver<-0
# }else{
#   if(length(DELX)==1){widthxvector<-matrix(DELX, 1, NX)
#   
#   }else{
#     widthxvector<-cbind(DELX)}
#   
#   distancetoriver<-vector(mode="numeric", length=NX)
#   distancetoriver[1]<- -0.5*widthxvector[1]
#   #    widthxvector[1]<- -widthxvector[1]*1/2
#   for(i in 2:NX){
#     distancetoriver[i]<-0.5*widthxvector[i]+0.5*widthxvector[i-1]+distancetoriver[i-1]}
# }
# distancetoriver[1]<-0
# 
# # define widthxvector here when no constant width is used
# widthyvector<-matrix(DELY, 1, NY) # define widthyvector here when no constant width is used

define_gwt_gridinput <- function (NX, NY, NRBL, widthxvector, widthyvector, NB=NULL, XB1=NULL, YB1=NULL, XB2=NULL, YB2=NULL, XB3=NULL, YB3=NULL, XB4=NULL, YB4=NULL) {
  firstline <- c(NX,NY, NRBL)
  fourthline <- NB
  rest <- cbind(c(XB1,XB2,XB3,XB4),c(YB1,YB2,YB3,YB4))
  write.table(t(firstline),"gwt_gridinput",row.names=FALSE,col.names=FALSE,sep=",")
  write.table(t(widthxvector),"gwt_gridinput",sep=",",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write.table(t(widthyvector),"gwt_gridinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write.table(c(fourthline),"gwt_gridinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write(paste(rest,sep=""),"gwt_gridinput",ncolumns=2,
        append=TRUE,sep=",")
}
#define_gwt_gridinput(NX, NY, NRBL, widthxvector, widthyvector)#, NB, XB1, YB1, XB2, YB2, XB3, YB3, XB4, YB4)

define_gwt_timestepinput <- function (ITIM, DELT, RECH, hriver) {
  firstline <- c(ITIM, DELT)
  write.table(t(firstline),"gwt_timestepinput", row.names=FALSE,col.names=FALSE,sep=",")
  write.table(format(t(RECH),digits=2),"gwt_timestepinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write.table(format(t(hriver),digits=2),"gwt_timestepinput", row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
}

#define_gwt_timestepinput(ITIM, DELT, RECH, hriver=hriver)

define_gwt_hydroinput <- function (headvector, bottomvector,ksat,spec_y) {
  firstline <- paste(headvector)
  secondline <- paste(bottomvector)
  thirdline <- paste(NX*NY,"*",ksat,sep="")
  fifthline <- paste(NX*NY,"*",spec_y,sep="")
  write.table(t(firstline),"gwt_hydroinput",row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE)
  write.table(t(secondline),"gwt_hydroinput",row.names=FALSE,col.names=FALSE,sep=",",
              append=TRUE,quote=FALSE)
  write.table(c(thirdline,thirdline,fifthline),"gwt_hydroinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  #  write.table(c(firstline,secondline,thirdline,thirdline,fifthline),"testgwt_hydroinput",row.names=FALSE,col.names=FALSE,
  #      append=TRUE,quote=FALSE)
  #  write(paste(RECH,"testgwt_gridinput"),ncolumns=NX*NY,
  #     append=TRUE,sep=",")
}


#ksat<-soilpar$K_s
#define_gwt_hydroinput(init_heads, bottom, "H Clay", ksat)                                                     

#define_gwt_gridinput(NX, NY, NRBL, widthxvector, widthyvector)#, NB, XB1, YB1, XB2, YB2, XB3, YB3, XB4, YB4)
#define_gwt_hydroinput(headvector, bottomvector, stype, ksat)                                                     
#define_gwt_timestepinput(ITIM, DELT, RECH)


define_gwt_riverinput <- function(NRBL,Ariver,criver,IXR,IXY) {
  rivermatrix <- matrix(0,NRBL,4) 
  rivermatrix[,3] <- ifelse(length(Ariver==1),rep(Ariver,NRBL),Ariver)
  rivermatrix[,4] <- ifelse(length(criver==1),rep(criver,NRBL),criver)
  #coordinates:
  rivermatrix[,1] <- IXR
  rivermatrix[,2] <- IXY
  write.table(rivermatrix,"gwt_riverinput", row.names=FALSE,col.names=FALSE,sep=",")
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


