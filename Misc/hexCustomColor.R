library(plotrix)
library(hexbin)

#creates a scale of colors
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

#generates data for three variables
dat=data.frame( x = c(rep(1:10,3)), y = c(rep(1:10,3)), z = c(rep(1:10,3)))

#generates hexbin with the x and y variables
hbin<-hexbin(dat$x, dat$y, xbins=10, IDs=TRUE)

#sums points falling inside bin
SumHexBin<-data.frame(sums=hexTapply(hbin, dat$z, sum))

#do color scale based on values of points in a third variable
cols <- myColorRamp(c("white","green","yellow", "red"), SumHexBin$sums)

## setup coordinate system of the plot
P <- plot(hbin, type="n",legend=FALSE)# asp=1

##add hexagons (in the proper viewport):
pushHexport(P$plot.vp)

#plots hexbins based on colors of third column
grid.hexagons(hbin, style= "lattice", border = gray(.9), pen = cols,  minarea = 1, maxarea = 1)
