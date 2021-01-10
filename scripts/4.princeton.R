rm(list = ls())
gc()
library(ncdf4)
library(plyr)
library(stringr)
library(lubridate)
library(PCICt)
library(sp)
library(raster)
library(fields)
library(abind)
library(rgdal)
library(ggplot2)


#1.princeton interp
dirnames <- dir(path = 'F:/princeton/prcp',full.names = TRUE)

output <- list()
for (i in seq_along(dirnames)){
  ncdata <- nc_open(filename = dirnames[i])
  lon <- ncvar_get(nc = ncdata ,varid = 'lon')
  lat <- ncvar_get(nc = ncdata , varid = 'lat')
  time <- ncvar_get(nc = ncdata, varid = 'time')
  time <- as.PCICt(x = '1948-01-01 00:00:0.0',cal = 'gregorian') + time*60
  prcp <- ncvar_get(nc = ncdata, varid = 'prcp')
  nc_close(ncdata)
  
  loc <- c(which(lon > 180),which(lon < 180))
  lon <- lon[loc]
  lon[lon > 180] <- lon[lon >180] - 360
  prcp <- prcp[loc,,]
  
  newdat <- list()
  for(j in 1:dim(prcp)[3]){
    newdat[[j]] <- interp.surface.grid(obj = list(x = lon ,y = lat ,z = prcp[,,j]),
                                       grid.list = list(x = seq(-179.75,179.75,0.5),
                                                        y = seq(-89.75,89.75,0.5)))$z
  }
  
  newdat <- abind(newdat,along = 3)
  attr(newdat , 'dimnames') <- NULL
  ym <- str_sub(string = time,start = 1,end = 7)
  monthly <- apply(X = newdat ,MARGIN = c(1,2),function(x){
    aggregate( x = x,by = list(ym),FUN = mean, na.rm = TRUE)$x
  })
  output[[i]] <- molnthly
}

prcp <- output
save(prcp,file = 'F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/month_prcp_Princeton.RData')


#change variable




#2.Z interp
dem <- readGDAL('F:/transfer/SyResearch/1.Attribution of global dryland areas/dem/dem_global_gtopo30.tif')
coord <- coordinates(dem)            # fetch the lonlat of tif
lon <- unique(coord[,1])             # get the unique
lat <- unique(coord[,2])
image(matrix(dem$band1, length(lon)))
z <- matrix(dem$band1, length(lon))

dem_z <- interp.surface.grid(obj = list(x = lon ,y = lat ,z = z ),
                             grid.list = list(x = seq(-179.75,179.75,0.5),
                                              y = seq(-89.75,89.75,0.5)))
# Z <- dem_z$z
save(Z, file = 'E:/SyResearch/1.Attribution of global dryland areas/output data/Z_0.5.RData')
lat <- dem_z$y
Mon <- rep(1:12, each = 69)
Z    <- array(data = rep(Z, times = 12), dim = c(720,360,12))
Z <- aperm(Z, c(3,1,2))
save(Z, file = 'F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/Z_i0.5.RData')



#3.cal pet of princeton

{
  load("F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/month_dswrf_Princeton.RData")
  load("F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/month_shum_Princeton.RData")
  load("F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/month_tmax_Princeton.RData")
  load("F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/month_tmin_Princeton.RData")
  load("F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/month_wind_Princeton.RData")
  load("F:/transfer/SyResearch/1.Attribution of global dryland areas/output data/Z_i0.5.RData")
}

Lonlat <- expand.grid(x = seq(-179.75,179.75,0.5),   y = seq(-89.75,89.75,0.5))
Lat <- array(rep(Lonlat$y, times = 12), dim = c(720, 360, 12))
Lat <- aperm(Lat, c(3, 1, 2))

Mon <- array(rep(1:12, each = 720*360), dim = c(720, 360, 12))
Mon <- aperm(Mon, c(3,1,2))


PET <- llply(1:length(dswrf), function(i){
  Penman_hist(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,
              Rs = dswrf[[i]]*24*60*60/1e6,
              WS = wind[[i]], SHum = shum[[i]],
              Mon = Mon, Z = Z, lat = Lat)
}, .progress = 'text')


PET_annual_wRs <- laply(PET,function(x){
  apply(x, c(2,3), mean, na.rm = T)
})
attr(PET_annual_wRs, 'dimnames') <- NULL
PET_annual_wRs <- aperm(PET_annual_wRs, c(2,3,1))

save(PET_annual_wRs, file = 'G:/1.dyland/output data/AI_his/preandpet/princeton/PET_annualmean.RData')


# PET_nRs<- llply(1:length(dswrf), function(i){
#   Penman_hist(Tmin =  tmin[[i]] - 273.16, Tmax = tmax[[i]]- 273.16,
#               WS = wind[[i]], SHum = shum[[i]],
#               Mon = Mon, Z = Z, lat = Lat)
# }, .progress = 'text')
# 
# 
# PET_annual_nRs <- laply(PET_nRs,function(x){
#   apply(x,c(2,3), sum, na.rm = T)
# })
# attr(PET_annual_nRs, 'dimnames') <- NULL
# 
# 
# save(PET_annual_nRs, file = 'E:/SyResearch/1.Attribution of global dryland areas/output data/PET_annual_nRs.RData')




#4.cal AI of princeton
{
  load("G:/1.dyland/output data/AI_his/preandpet/princeton/prcp_annualmean.RData")
  load("G:/1.dyland/output data/AI_his/preandpet/princeton/PET_annualmean.RData")
}

# f <- function(x)x*24*3600
# prcp1 <- lapply(prcp,f)
# prcp_annual <- lapply(prcp1, function(x){
#   apply(x, c(2,3), sum, na.rm = T)
# })
# PET_annual_wRs <- aperm(PET_annual_wRs, c(2,3,1))
# prcp_annual <-  abind(prcp_annual,along = 3)
# attr(prcp_annual,"dimnames") <- NULL
# prcp_annual <- aperm(prcp_annual,c(1,2,3))

prcp_mean <- prcp_mean*24*3600
# PET_annual_wRs[PET_annual_wRs < 0] <- NA
load("G:/1.dyland/output data/AI_his/1958-2014/AI_cru57.RData")
prcp_mean_ns <- laply(1:69, function(i){
  x <- prcp_mean[,,i]
  x[is.na(AI_cru57[,,1])] <- NA
  x
})
attr(prcp_mean_ns, 'dimnames') <- NULL
prcp_mean_ns <- aperm(prcp_mean_ns, c(2,3,1))


PET_annual_wRs <- laply(1:69, function(i){
  x1 <- PET_annual_wRs[,,i]
  x1[is.na(AI_cru57[,,1])] <- NA
  x1
})
attr(PET_annual_wRs, 'dimnames') <- NULL
PET_annual_wRs <- aperm(PET_annual_wRs, c(2,3,1))

AI_princeton <- prcp_mean_ns/PET_annual_wRs

save(PET_annual_wRs, file = "G:/1.dyland/output data/AI_his/preandpet/princeton/PET_annualmeannosea.RData")
save(prcp_mean_ns,file = "G:/1.dyland/output data/AI_his/preandpet/princeton/prcp_mean_ns.RData")
save(AI_princeton, file = "G:/1.dyland/output data/AI_his/all/AI_princeton1948-2016.RData")


AI_princeton57 <- AI_princeton[,,11:67]
save(AI_princeton57,file = "G:/1.dyland/output data/AI_his/1958-2014/AI_princeton57.RData")


AI_princeton30 <- AI_princeton[,,14:43]
save(AI_princeton30, file = "G:/1.dyland/output data/AI_his/1961-1990/AI_princeton30.RData")

# 
# #5.delete sea
# #5.1 delete sea
# cru <- nc_open("G:/Observation/CRU/cru_ts4.04.1901.2019.pet.dat.nc")
# cru_pet <- ncvar_get(cru, varid = 'pet')
# load("E:/Research/output/AI_Princeton/AI_princeton_wRs.RData")
# 
# 
# hist(cru_pet) # Note that the unit of Cru_pet is mm/day
# 
# library(plyr)
# AI2 <- laply(1:69, function(i){
#   x2 <- AI_princeton_wRs[,,i]
#   x2[is.na(cru_pet[,,1])] <- NA
#   x2
# })
# attr(AI2, 'dimnames') <- NULL
# summary(AI2)
# AI2[AI2 > 2] <- 2
# AI2[AI2 < 0] <- 0
# image(AI2[1,,],)
# hist(AI2)


#AI1 <- apply(AI_princeton_wRs, c(1,2), function(x) x[is.na(cru_pet[,,1])] <- NA)
#this method not get success

#6.cal area of princeton
load("G:/1.dyland/output data/AI_his/1958-2014/AI_princeton57.RData")
lon = seq(-179.75, 179.75, 0.5)
lat = seq(-89.75, 89.75, 0.5)
lonlat <- expand.grid(lon,lat)
Area <- F_area(lonlat)
princeton_area57 <- matrix(NA,57)
for(i in 1:57){
  princeton_area57[i] <- sum(Area[AI_princeton57[,,i] <0.65],na.rm = T)
}

save(princeton_area57, file = "G:/1.dyland/output data/area_his/1958-2014/princeton_area57.RData")



load("G:/1.dyland/output data/AI_his/1961-1990/AI_princeton30.RData")
lon = seq(-179.75, 179.75, 0.5)
lat = seq(-89.75, 89.75, 0.5)
lonlat <- expand.grid(lon,lat)
Area <- F_area(lonlat)
princeton_area30 <- matrix(NA,30)
for(i in 1:30){
  princeton_area30[i] <- sum(Area[AI_princeton30[,,i] <0.65],na.rm = T)
}

save(princeton_area30, file = "G:/1.dyland/output data/area_his/1961-1990/princeton_area30.RData")






#5.2 trend figure
load("G:/1.dyland/output data/area_his/1958-2014/princeton_area57.RData")


Prin <- data.frame(Years = 1958:2014,Area = princeton_area57)
p1<-ggplot(Prin,aes(x = Years,y= Area))+geom_point(shape=17,size=2.4,colour = 'red')+
  stat_smooth(method = lm,level = 0.99,size = 1.2)+scale_color_manual(values = '#CC6666')
p1+theme_bw()+ggtitle("(a)")



#================================================================================end

lon <- ncvar_get(cru,varid = 'lon')
lat <- ncvar_get(cru,varid = 'lat')
lonlat <- expand.grid(lon,lat)

Area <- F_area(lonlat)
Areas <- matrix(NA,69)
for(i in 1:69){
  Areas[i] <- sum(Area[AI2[i,,] <0.65],na.rm = T)
}
Prin <- data.frame(Years = 1948:2016,Area = Areas)
p1<-ggplot(Prin,aes(x = Years,y= Area))+geom_point(shape=17,size=2.4,colour = 'red')+
  stat_smooth(method = lm,level = 0.99,size = 1.2)+scale_color_manual(values = '#CC6666')
p1+theme_bw()+ggtitle("(a)")+mytheme

library(trend)
sens.slope(Areas)


AI58 <- AI2[11:68,,]
Area <- F_area(lonlat)
Areas <- matrix(NA,58)
for(i in 1:58){
  Areas[i] <- sum(Area[AI58[i,,] <0.65],na.rm = T)
}
Prin <- data.frame(Years = 1958:2015,Area = Areas)
p1<-ggplot(Prin,aes(x = Years,y= Area))+geom_point(shape=17,size=2.4,colour = 'red')+
  stat_smooth(method = lm,level = 0.99,size = 1.2)+scale_color_manual(values = '#CC6666')
p1+theme_bw()+ggtitle("(a)")+mytheme



#3. no Rs
load("E:/Research/output/AI_Princeton/AI_princeton_nRs.RData")
AI3 <- laply(1:69, function(i){
  x3 <- AI_princeton_nRs[,,i]
  x3[is.na(cru_pet[,,1])] <- NA
  x3
})
attr(AI3, 'dimnames') <- NULL
summary(AI3)
AI3[AI3 > 2] <- 2
AI3[AI3 < 0] <- 0

for(i in 1:69){
  Areas[i] <- sum(Area[AI3[i,,] <0.65],na.rm = T)
}
Prin <- data.frame(Years = 1948:2016,Area = Areas)
p1<-ggplot(Prin,aes(x = Years,y= Area))+geom_point(shape=17,size=2.4,colour = 'red')+
  stat_smooth(method = lm,level = 0.99,size = 1.2)+scale_color_manual(values = '#CC6666')
p1+theme_bw()+ggtitle("(a)")+mytheme

library(trend)
sens.slope(Areas)

AI3_58 <- AI3[11:68,,]
Area <- F_area(lonlat)
Areas <- matrix(NA,58)
for(i in 1:58){
  Areas[i] <- sum(Area[AI3_58[i,,] <0.65],na.rm = T)
}
Prin <- data.frame(Years = 1958:2015,Area = Areas)
p1<-ggplot(Prin,aes(x = Years,y= Area))+geom_point(shape=17,size=2.4,colour = 'red')+
  stat_smooth(method = lm,level = 0.99,size = 1.2)+scale_color_manual(values = '#CC6666')
p1+theme_bw()+ggtitle("(a)")+mytheme





























