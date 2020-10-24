library(ncdf4)
sample <- nc_open('E:/SyResearch/1.Attribution of global dryland areas/sample/cru_ts4.04.1901.2019.pet.dat.nc')
sample <- ncvar_get(sample,varid = 'pet')
image(sample[,,565])
cor.test(PET2[1,,],sample[,,565])
hist(PET2[1,,])
hist(sample[,,556])
sim1<- PET[1,,]
sim2<- PET2[1,,]
obs <- sample[,,565]

# image(sim)
# image(obs)
# cor.test(sim, obs)
plot(sim, obs)
plot(sim2, obs)
cor.test(sim2, obs)
