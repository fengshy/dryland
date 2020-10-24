setwd('E:/SyResearch/1.Attribution of global dryland areas/codes/')
rm(list = ls())
gc()
library(ncdf4)
library(plyr)
library(stringr)
library(lubridate)
library(PCICt)
library(sp)
library(fields)
library(abind)

dirnames <- dir(path = 'G:/princeton/tmin',full.names = TRUE)
#ncdata <- nc_open(filename = dirnames[72])
#names(ncdata$var)
#names(ncdata$dim)

#lon <- ncvar_get(nc = ncdata ,varid = 'lon')
#lat <- ncvar_get(nc = ncdata , varid = 'lat')
#time <- ncvar_get(nc = ncdata, varid = 'time')

output <- list()
for (i in seq_along(dirnames)){
  ncdata <- nc_open(filename = dirnames[i])
  lon <- ncvar_get(nc = ncdata ,varid = 'lon')
  lat <- ncvar_get(nc = ncdata , varid = 'lat')
  time <- ncvar_get(nc = ncdata, varid = 'time')
  time <- as.PCICt(x = '1948-01-01 00:00:0.0',cal = 'gregorian') + time*60
  tmin <- ncvar_get(nc = ncdata, varid = 'tmin')
  nc_close(ncdata)

  loc <- c(which(lon > 180),which(lon < 180))
  lon <- lon[loc]
  lon[lon > 180] <- lon[lon >180] - 360
  tmin <- tmin[loc,,]

  newdat <- list()
  for(j in 1:dim(tmin)[3]){
    newdat[[j]] <- interp.surface.grid(obj = list(x = lon ,y = lat ,z = tmin[,,j]),
                                       grid.list = list(x = seq(-179.75,179.75,0.5),
                                                        y = seq(-89.75,89.75,0.5)))$z
    print(j)
  }

  newdat <- abind(newdat,along = 3)
  attr(newdat , 'dimnames') <- NULL
  ym <- str_sub(string = time,start = 1,end = 7)
  monthly <- apply(X = newdat ,MARGIN = c(1,2),function(x){
    aggregate( x = x,by = list(ym),FUN = mean, na.rm = TRUE)$x
  })
  output[[i]] <- monthly
}

tmin <- output
save(tmin,file = 'E:/SyResearch/1.Attribution of global dryland areas/output data/month_tmin_Princeton.RData')









