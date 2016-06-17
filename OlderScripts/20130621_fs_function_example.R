# testing the new fs function example

setwd("x:/vervoort/research/ecohydrology/2dmodelling")
rdir <- "x:/vervoort/research/rcode/ecohydrology/2dmodelling"


source(paste(rdir,"20120724_FluxfunctionsforElise.R",sep="/"))
source(paste(rdir,"soilfunction.R",sep="/"))
source(paste(rdir,"vegfunction.R",sep="/"))

soilpar <- Soil("L Med Clay")
vegpar <- Veg(vtype="TreesDR", soilpar=soilpar)
DR <- vegpar$DR
# This is a key variable based on Vervoort and van der Zee (2012)
fs <- seq(0.2,0.7,length=5)
fs <- c(0.25,0.5)
Z <- seq(700,300)
Zmean <- 500


plot(RWU(Z,Zmean,fs[1]),Z,type="l",xlab="Root water uptake function",ylab="depth (cm)",
      xlim=c(0,8),lwd=2,cex.axis=1.2,cex.lab=1.2,font.lab=2)
for (i in 2:length(fs)) {
  lines(RWU(Z,Zmean,fs[i]),Z,lty=i,col=i,lwd=2)
}

# this shows that the new root function already includes the effect of anoxia
# if water goes above a certain Z value the root water uptake decreases
plot(Rc_B(Z,vegpar$c1,vegpar$Zr,Zmean,fs[1]),Z,type="l",xlab="Root water uptake function",ylab="depth (cm)",
      xlim=c(0,8),lwd=2,cex.axis=1.2,cex.lab=1.2,font.lab=2)
for (i in 2:length(fs)) {
  lines(Rc_B(Z,vegpar$c1,vegpar$Zr,Zmean,fs[i]),Z,lty=i,col=i,lwd=2)
}

vegpar$c1 <- 1.5
plot(U(z1=Z,z2=vegpar$Zr,c1=vegpar$c1),Z)



legend("topright",paste("fs =",fs),lty=1:5,col=1:5,lwd=2,cex=1.2)