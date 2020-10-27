rm(list = ls())
gc()

library(ncdf4)
library(sp)
library(rgdal)
library(plyr)
library(fields)

{
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_dswrf_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_shum_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_tmax_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_tmin_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_wind_Princeton.RData")
  load("E:/SyResearch/1.Attribution of global dryland areas/output data/Z_i0.5.RData")
}

Lonlat <- expand.grid(x = seq(-179.75,179.75,0.5),   y = seq(-89.75,89.75,0.5))
Lat <- array(rep(Lonlat$y, times = 12), dim = c(720, 360, 12))
Lat <- aperm(Lat, c(3, 1, 2))

Mon <- array(rep(1:12, each = 720*360), dim = c(720, 360, 12))
Mon <- aperm(Mon, c(3,1,2))


PET<- llply(1:length(dswrf), function(i){
  Penman_hist(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,
                        Rs = dswrf[[i]]*24*60*60/1e6,
                        WS = wind[[i]], SHum = shum[[i]],
                        Mon = Mon, Z = Z, lat = Lat)
}, .progress = 'text')


PET_annual_wRs <- laply(PET,function(x){
  apply(x,c(2,3), sum, na.rm = T)
})
attr(PET_annual_wRs, 'dimnames') <- NULL


PET_nRs<- llply(1:length(dswrf), function(i){
  Penman_hist(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,
              WS = wind[[i]], SHum = shum[[i]],
              Mon = Mon, Z = Z, lat = Lat)
}, .progress = 'text')


PET_annual_nRs <- laply(PET_nRs,function(x){
  apply(x,c(2,3), sum, na.rm = T)
})
attr(PET_annual_nRs, 'dimnames') <- NULL

save(PET_annual_wRs, file = 'E:/SyResearch/1.Attribution of global dryland areas/output data/PET_annual_wRs.RData')
save(PET_annual_nRs, file = 'E:/SyResearch/1.Attribution of global dryland areas/output data/PET_annual_nRs.RData')
