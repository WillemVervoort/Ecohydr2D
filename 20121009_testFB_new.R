# test FB_new

today <- format(Sys.Date(),"%Y%m%d")
#setwd("C:/Documents and Settings/Administrator/Desktop/Hons/3-My Work/4 - R files")
setwd("x:/research/ecohydrology/2dmodelling")
#rdir <- "C:/Documents and Settings/Administrator/Desktop/Hons/3-My Work/4 - R files"
rdir <- "x:/research/rcode/ecohydrology/2dmodelling"


source(paste(rdir,"20120724_FluxfunctionsforElise.R",sep="/"))


s <- seq(0.2,1,length=100)

soilpar <- Soil("L Med Clay")
vegpar <- Veg(vtype="TreesDR", soilpar=soilpar)
DR <- vegpar$DR
# This is a key variable based on Vervoort and van der Zee (2012)
vegpar$fs <- 0.3

Z <- 1075
Zmean <- 1075

# run FB_new
foo <- data.frame(matrix(ncol=5,nrow=length(s)))

for (i in 1:length(s)) {
  foo[i,] <- FB_new(s[i],soilpar,vegpar, Z, Zmean,Z.prev=1080)
}

head(foo)
