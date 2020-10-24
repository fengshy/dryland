setwd('E:/SyResearch/1.Attribution of global dryland areas/codes/')
rm(list = ls())
gc()
library(ncdf4)
library(sp)
library(rgdal)
library(plyr)
library(fields)
dem <- readGDAL('E:/SyResearch/1.Attribution of global dryland areas/dem/dem_global_gtopo30.tif')
coord <- coordinates(dem)            # fetch the lonlat of tif
lon <- unique(coord[,1])             # get the unique
lat <- unique(coord[,2])
rm(coord)
#image(matrix(dem$band1, length(lon)))
z <- matrix(dem$band1, length(lon))

dem_z <- interp.surface.grid(obj = list(x = lon ,y = lat ,z = z ),
                                     grid.list = list(x = seq(-179.75,179.75,0.5),
                                                      y = seq(-89.75,89.75,0.5)))

Z <- dem_z$z
image(Z)
save(Z, file = 'E:/SyResearch/1.Attribution of global dryland areas/output data/Z_0.5.RData')
lat <- dem_z$y
Mon <- rep(1:12, each = 69)
