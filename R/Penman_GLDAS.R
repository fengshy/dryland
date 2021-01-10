# Penman-Monteith Evapotranspiration (FAO-56 Method), 2020.10.27
# Comparing with the example in FAO handbook, results are right.
# Comparing with the results based on the package of evapotranspiration, results are right.

# Two dataset (Princeton and CPC) have Rs and specific humidity (SHum), but no RH and Rn.

#' saturation vapor pressure
#' @param T temperature, in degC
#' @return es, in kPa
#' @export
cal_es <- function(T) {
  0.6108*exp(17.27*T/(T+237.3))
}

#' @export
Penman_his <- function(Tm, WS, SHum, Rns, Rnl, Z){
  # 1. Calculate Rn
  Rn = Rns + Rnl

  #2. Calculate es and ea
  press = 101.325
  es <- cal_es(Tm)
  ea <- SHum * press / (0.378 * SHum + 0.622)


  # 3. Calculate PET
  delta <- 4098*cal_es(Tm)/(Tm+237.3)^2
  gama <- 0.000665*101.3*((293-0.0065*Z)/293)^5.26
  PET <- (0.408*delta*Rn+gama*900/(Tm+273)*WS*(es-ea))/(delta+gama*(1+0.34*WS))
  return(PET)
}
