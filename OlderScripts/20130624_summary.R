# read and plot summary

setwd("x:/vervoort/research/ecohydrology/2dmodelling")

summar <- read.csv("20130624_summary.csv",skip=1)

summar_TG <- stack(summar[,c(1,3,5,7,9)])
summar_stress <- stack(summar[,c(2,4,6,8,10)])

summar1 <- data.frame(grid = rep(seq(5,95,length=10),5), TG = summar_TG,
			stress = summar_stress, id = c(rep("Base",10),rep("fs = 0.25",10),
			rep("No river",10), rep("No rain",10), rep("Out of phase",10)))

require(lattice)
barchart(TG.values ~ as.factor(grid) | id, data = summar1, horizontal=FALSE,
	ylab=list("Total Groundwater Transpiration",font=2,cex=1.2),
    scales=list(font=2,cex=1.2),xlab=list("Distance (m)",font=2,cex=1.2),as.Table=T)

barchart(stress.values ~ as.factor(grid) | id, data = summar1, horizontal=FALSE,
	ylab=list("Average Dynamic Stress",font=2,cex=1.2),
    scales=list(font=2,cex=1.2),xlab=list("Distance (m)",font=2,cex=1.2),
ylim=c(0,1),as.Table=T)