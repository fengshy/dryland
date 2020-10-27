setwd('E:\\Research\\Attribution of global dryland areas\\codes')
rm(list = ls())
gc()
library(ncdf4)
library(plyr)
library(stringr)
library(lubridate)
library(PCICt)
library(sp)

dirnames <- dir(path = 'E:\\princeton\\data', full.names = TRUE)
ncdata <- nc_open(filename = dirnames[1])
names(ncdata$dim)
names(ncdata$var)
ncdata$var$dlwrf$units
ncdata$dim$time$units

lon <- ncvar_get(nc = ncdata, varid = 'lon')
lat <- ncvar_get(nc = ncdata, varid = 'lat')
time <- ncvar_get(nc = ncdata, varid = 'time')

for (i in 1:length(time)){
  as.PCICt(x = '1948-01-01 00:00:0.0', cal = 'gregorian') + time[i]*60
}


time <- as.PCICt(x = '1948-01-01 00:00:0.0', cal = 'gregorian') + time*60
time <- format.Date(x = time, format = '%Y-%m-%d')
dlwrf <- ncvar_get(nc = ncdata, varid = 'dlwrf')
nc_close(ncdata)
grid1 <- SpatialPixelsDataFrame(points = expand.grid(lon, lat),
                                data = data.frame(val = c(dlwrf[,,1])),
                                proj4string = CRS('+proj=longlat +ellps=WGS84'))

loc <- c(which(lon > 180), which(lon < 180))
lon.new <- lon[loc]
lon.new[lon.new > 180] <- lon.new[lon.new > 180] - 360

dlwrf.new <- dlwrf[loc, , ]

grid2 <- SpatialPixelsDataFrame(points = expand.grid(lon.new, lat),
                                data = data.frame(val = c(dlwrf.new[,,1])),
                                proj4string = CRS('+proj=longlat +ellps=WGS84'))

grid3 <- SpatialPixelsDataFrame(points = expand.grid(lon.new, lat),
                                data = data.frame(val = c(dlwrf[,,1])),
                                proj4string = CRS('+proj=longlat +ellps=WGS84'))

spplot(grid3)
