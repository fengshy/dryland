library(raster)
library(magrittr)
library(rgdal)

library(nctools)
library(matrixStats)
library(ncdf4)
library(lubridate)


# "minutes since 1948-01-01 00:00"
nc_day2month <- function(file, outfile) {
    fid = nc_open(file)
    dates = ymd("1948-01-01") + fid$dim$time$vals/1440

    b <- brick(file)
    arr <- as.array(b)
    mat = apply_3d(arr, 3, rowMeans2, by = month(dates))

    range = b@extent # S4 class '@', S3 $, spatial data
    b2 = brick(mat, range@xmin, range@xmax, range@ymin, range@ymax)
    writeRaster(b2, outfile)
}

files <- dir("e:/princeton/", "*.nc", full.names = TRUE)
file <- "e:/princeton/dlwrf_daily_1948-1948.nc"
nc_day2month(file, "output_1948.nc")

## version 1
for(file in files) {
  print(file)
}

## version 2
for(i in seq_along(files)) {
  file = files[i]
  outfile = paste('dlwrf','_mon_',1948+i-1,".nc",sep = '')

  print(file)
  nc_day2month(file, outfile)
  # if(i=1,break)
  # print(file)
}


# b3 = raster("output.nc")
# fid <- nc_open("output.nc")
# x = ncvar_get(fid, "variable")
