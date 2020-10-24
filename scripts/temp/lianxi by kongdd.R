library(raster)
library(magrittr)
library(rgdal)

library(nctools)
library(matrixStats)
library(ncdf4)
library(lubridate)


setwd('E:/princeton/')
# nc_date(file)
daily2mon <- function(i,variable){

"minutes since 1948-01-01 00:00"
  file <- paste("E:/princeton/",variable,'_daily_',1948+i-1,'-',1948+i-1,".nc",sep = "")
  fid = nc_open(file)
  h <- paste(1948+i-1,'-01-01')
  dates = ymd(h) + fid$dim$time$vals/1440
  
  b <- brick(file)
  arr <- as.array(b)
  
  mat = apply_3d(arr, 3, rowMeans2, by = month(dates))
  
  range = b@extent
  b2 = brick(mat, range@xmin, range@xmax, range@ymin, range@ymax)
  k <- paste(variable,'_mon_',1948+i-1,'.nc')
  writeRaster(b2, 'k')
}
"minutes since 1948-01-01 00:00"
out1950 <- daily2mon(i=3,variable ='dlwrf' )







b3 = raster("output.nc")
# fid <- nc_open("output.nc")
# x = ncvar_get(fid, "variable")
