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

{
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_dswrf_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_shum_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_tmax_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_tmin_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_wind_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/Z_0.5.RData")
}
Z    <- array(data = rep(Z, times = 12), dim = c(720,360,12))
Z <- aperm(Z, c(3,1,2))

Lonlat <- expand.grid(x = seq(-179.75,179.75,0.5),   y = seq(-89.75,89.75,0.5))
Lat <- array(rep(Lonlat$y, times = 12), dim = c(720, 360, 12))
Lat <- aperm(Lat, c(3, 1, 2))

Mon <- array(rep(1:12, each = 720*360), dim = c(720, 360, 12))
Mon <- aperm(Mon, c(3,1,2))

#i <- 1
# PET <- Penman_history(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,  WS = wind[[i]], SHum = shum[[i]],
#                       Rs = dswrf[[i]], Mon = Mon, Z = Z, lat = Lat)
# PET <- laply(1, function(i){
  PET <- Penman_history(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,
                        Rs = dswrf[[i]]*24*60*60/1e6,
                        WS = wind[[i]], SHum = shum[[i]],
                        Mon = Mon, Z = Z, lat = Lat)
  {
    x = PET[1,,]
    #x[x < 0] = NA
    levelplot(x*365)
    #levelplot(dswrf[[i]][1,,]*365*24*60*60/1e6)
  }
  # PET
# }, .progress = 'text')
# attr(PET,'dimnames') <- NULL
