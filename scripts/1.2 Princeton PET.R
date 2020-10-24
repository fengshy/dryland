setwd('E:/SyResearch/1.Attribution of global dryland areas/codes/')
rm(list = ls())
gc()


load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_dswrf_Princeton.RData")
load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_shum_Princeton.RData")
load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_tmax_Princeton.RData")
load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_tmin_Princeton.RData")
load("E:/SyResearch/1.Attribution of global dryland areas/output data/month_wind_Princeton.RData")
load("E:/SyResearch/1.Attribution of global dryland areas/output data/Z_0.5.RData")

Z    <- array(data = rep(Z, times = 12), dim = c(720,360,12))
Z <- aperm(Z, c(3,1,2))

Lonlat <- expand.grid(x = seq(-179.75,179.75,0.5),   y = seq(-89.75,89.75,0.5))
Lat <- array(rep(Lonlat$y, times = 12), dim = c(720, 360, 12))
Lat <- aperm(Lat, c(3, 1, 2))

Mon <- array(rep(1:12, each = 720*360), dim = c(720, 360, 12))
Mon <- aperm(Mon, c(3,1,2))

source('PM_feng.R')

#i <- 1
# PET <- Penman_history(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,  WS = wind[[i]], SHum = shum[[i]],
#                       Rs = dswrf[[i]], Mon = Mon, Z = Z, lat = Lat)
library(plyr)
PET <- laply(1:69, function(i){
  
  PET <- Penman_history_2(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,  WS = wind[[i]], SHum = shum[[i]],
                           Mon = Mon, Z = Z, lat = Lat)
  PET
}, .progress = 'text')



attr(PET,'dimnames') <- NULL
