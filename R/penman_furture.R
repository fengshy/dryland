
#use for CMIP6
#' @export
Penman_future <- function(Tmax,Tmin,WS,RH, SH, LH, Z){
  Tm <- (Tmax+Tmin)/2
  Rn <- SH + LH  # sensible heat and Latent heat flux 
  delta <- 4098*(0.6108*exp(17.27*Tm/(Tm+237.3)))/(Tm+237.3)^2
  gama <- 0.000665*101.3*((293-0.0065*Z)/293)^5.26
  es <- (0.6108*exp(17.27*Tmax/(Tmax+237.3))+0.6108*exp(17.27*Tmin/(Tmin+237.3)))/2
  ea <- RH/100*(0.6108*exp(17.27*Tmax/(Tmax+237.3))+0.6108*exp(17.27*Tmin/(Tmin+237.3)))/2
  PET <- (0.408*delta*Rn+gama*900/(Tm+273)*WS*(es-ea))/(delta+gama*(1+0.34*WS))
  return(PET)
}
# save(PET,file = 'E:/SyResearch/1.Attribution of global dryland areas/sample/PET_outofRs.R')

