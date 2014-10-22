library(irisws)
library(reshape2)
library(plyr)

set.seed(111)

n <- 5
ns <- seq.int(1,n)
stations <- list(name = letters[ns],
  loc=list(lat = 32 + rnorm(n), lon = -116 + rnorm(n)))

Dist <- combn(rev(ns), n-1, tabulate, nbins = n)

all.distances <- function(lats, lons){
  
  cbind(expand.grid(sta.Lat=lats, sta.Lon=lons), degreelen(lats)[,2:3])
}

tail(with(stations$loc, all.distances(lat,lon)))
