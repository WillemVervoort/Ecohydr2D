# This set of functions writes the input files for GWT.exe
# it assumes that GWT.exe is in a sub folder called /GWTmod

#' Generate the groundwater grid input data for GWT.exe
#'
#' \code{define_gwt_gridinput} generates the grid input data for GWT.exe
#'
#' This is a utility function specifically designed for GWT.exe It is mostly
#' used internally and therefore not directly used. Refer to GWT.exe
#' documentation. It assumes that GWT.exe is in a sub folder called /GWTmod
#'
#' @param NX Number of gridcells in the X direction
#' @param NY Number of gridcells in the y direction
#' @param NRBL Number of river blocks
#' @param widthxvector vector of widths (\code{length(NX)}) of the blocks in the x-direction
#' @param widthyvector vector of widths (\code{length(NY)}) of the blocks in the y-direction
#' @param NB (optional) number of blocks on the boundary of a polygon
#' @param XB1 (optional) X coordinate of the first point on the boundary of a polygon
#' @param YB1 (optional) Y coordinate of the first point on the boundary of a polygon
#' @param XB2 (optional) X coordinate of the second point on the boundary of a polygon
#' @param YB2 (optional) Y coordinate of the second point on the boundary of a polygon
#' @param XB3 (optional) X coordinate of the third point on the boundary of a polygon
#' @param YB3 (optional) Y coordinate of the third point on the boundary of a polygon
#' @param XB4 (optional) X coordinate of the fourth point on the boundary of a polygon
#' @param YB4 (optional) Y coordinate of the fourth point on the boundary of a polygon
#' @examples
#' define_gwt_gridinput(10,1,2,rep(10,12),10)
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
## test
# define_gwt_gridinput(10, 1, 2, rep(10,12), 10)

#'Generate the timestep input data for GWT.exe
#'
#'\code{define_gwt_timestepinput} generates the timestep input data for GWT.exe
#'
#'This is a utility function specifically designed for GWT.exe. It is mostly
#'used internally and therefore not directly used. Refer to GWT.exe
#'documentation. It assumes that GWT.exe is in a sub folder called /GWTmod
#'
#'@param ITIM timestep number for the groundwater model
#'@param DELT timestep size
#'@param RECH vector of recharge for each of the grid blocks \code{length(NX)}
#'@param hriver vector with height of water in the river blocks
#'  \code{length(NRBL)}
#'@examples
#'NX <- 10
#'NRBL <- 2
#'define_gwt_timestepinput(1,1,rep(0.1,NX),rep(2,NRBL))
define_gwt_timestepinput <- function (ITIM, DELT, RECH, hriver) {
  firstline <- c(ITIM, DELT)
  write.table(t(firstline),"GWTmod/gwt_timestepinput", row.names=FALSE,col.names=FALSE,sep=",")
  write.table(format(t(RECH),digits=2),"GWTmod/gwt_timestepinput",row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
  write.table(format(t(hriver),digits=2),"GWTmod/gwt_timestepinput", row.names=FALSE,col.names=FALSE,
              append=TRUE,quote=FALSE)
}


#' Generate the hydraulic data for GWT.exe
#'
#' \code{define_gwt_hydroinput} generates the hydraulic data for GWT.exe
#'
#' This is a utility function specifically designed for GWT.exe It is mostly
#' used internally and therefore not directly used. Refer to GWT.exe
#' documentation. It assumes that GWT.exe is in a sub folder called /GWTmod
#'
#' @param NX Number of gridcells in the X direction
#' @param NY Number of gridcells in the y direction
#' @param headvector vector of the groundwater head in each cell (\code{NX*NY})
#' @param bottomvector vector of the height of the aquifer bottom (\code{NX*NY})
#' @param ksat in the x and y direction
#' @param spec_y specific yield of the aquifer
#' @examples
#' NX <- 10
#' NY <- 1
#' define_gwt_hydroinput(NX,NY, rep(10,NX*NY), rep(-25,NX*NY),1.5,0.05)
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


#' Generate the river input data for GWT.exe
#'
#' \code{define_gwt_riverinput} generates the river input data for GWT.exe
#'
#' This is a utility function specifically designed for GWT.exe It is mostly
#' used internally and therefore not directly used. Refer to GWT.exe
#' documentation. It assumes that GWT.exe is in a sub folder called /GWTmod
#'
#' @param NRBL Number of river grid cells
#' @param Ariver Area of the river blocks
#' @param criver resistance between river and aquifer block
#' @param IXR x coordinates of the river blocks
#' @param IYR y coordinates of the river blocks
#' @examples
#' define_gwt_riverinput(2,200,0.01,c(1,12),c(1,1))
define_gwt_riverinput <- function(NRBL,Ariver,criver,IXR,IYR) {
  rivermatrix <- matrix(0,NRBL,4)
  rivermatrix[,3] <- Ariver
  rivermatrix[,4] <- criver
  #coordinates:
  rivermatrix[,1] <- IXR
  rivermatrix[,2] <- IYR
  write.table(rivermatrix,"GWTmod/gwt_riverinput", row.names=FALSE,col.names=FALSE,sep=",")
}
