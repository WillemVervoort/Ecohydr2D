# This set of functions writes the input files for GWT.exe
# it assumes that GWT.exe is in a sub folder called /GWTmod

# gwt grid input
define_gwt_gridinput <- function (NX, NY, NRBL, 
                                  widthxvector, widthyvector, NB=NULL, XB1=NULL,
                                  YB1=NULL, XB2=NULL, YB2=NULL, XB3=NULL, 
                                  YB3=NULL, XB4=NULL, YB4=NULL) {
  firstline <- c(NX,NY, NRBL)
  fourthline <- NB
  rest <- cbind(c(XB1,XB2,XB3,XB4),c(YB1,YB2,YB3,YB4))
  write.table(t(firstline),"GWTmod/gwt_gridinput",row.names=FALSE,col.names=FALSE,sep=",")
  write.table(t(widthxvector),"GWTmod/gwt_gridinput",sep=",",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write.table(t(widthyvector),"GWTmod/gwt_gridinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write.table(c(fourthline),"GWTmod/gwt_gridinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write(paste(rest,sep=""),"GWTmod/gwt_gridinput",ncolumns=2,
        append=TRUE,sep=",")
}
#define_gwt_gridinput(NX, NY, NRBL, widthxvector, widthyvector)#, NB, XB1, YB1, XB2, YB2, XB3, YB3, XB4, YB4)

# define gwt timestepinput file
define_gwt_timestepinput <- function (ITIM, DELT, RECH, hriver) {
  firstline <- c(ITIM, DELT)
  write.table(t(firstline),"GWTmod/gwt_timestepinput", row.names=FALSE,col.names=FALSE,sep=",")
  write.table(format(t(RECH),digits=2),"GWTmod/gwt_timestepinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write.table(format(t(hriver),digits=2),"GWTmod/gwt_timestepinput", row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
}

#define_gwt_timestepinput(ITIM, DELT, RECH, hriver=hriver)

define_gwt_hydroinput <- function (NX,NY,headvector, bottomvector,
                                   ksat,spec_y) {
  firstline <- paste(headvector)
  secondline <- paste(bottomvector)
  thirdline <- paste(NX*NY,"*",ksat,sep="")
  fifthline <- paste(NX*NY,"*",spec_y,sep="")
  write.table(t(firstline),"GWTmod/gwt_hydroinput",row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE)
  write.table(t(secondline),"GWTmod/gwt_hydroinput",row.names=FALSE,col.names=FALSE,sep=",",
              append=TRUE,quote=FALSE)
  write.table(c(thirdline,thirdline,fifthline),"GWTmod/gwt_hydroinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
}


# River input file
define_gwt_riverinput <- function(NRBL,Ariver,criver,IXR,IYR) {
  rivermatrix <- matrix(0,NRBL,4) 
  rivermatrix[,3] <- Ariver
  rivermatrix[,4] <- criver
  #coordinates:
  rivermatrix[,1] <- IXR
  rivermatrix[,2] <- IYR
  write.table(rivermatrix,"GWTmod/gwt_riverinput", row.names=FALSE,col.names=FALSE,sep=",")
}
