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
Penman_hist <- function(Tmin, Tmax, WS, SHum, Rs = NULL, Mon, Z, lat){
  Tm <- (Tmin + Tmax)/2
  # 1. Calculate Ra
  j <- (Mon-1)*30+ 15
  #j <- i if cal daily,delete line19 and take line 20 in function(),delete Mon
  Gsc <- 0.0820
  dr <- 1+0.033*cos(2*pi/365*j)
  x <- pi/180*lat
  y <- 0.409*sin(2*pi/365*j-1.39)

  temp = -tan(x) * tan(y)
  temp[temp < -1] <- -1
  temp[temp >  1] <- 1
  ws <- acos(temp)

  Ra <- 24*60/pi*Gsc*dr*(ws*sin(x)*sin(y)+cos(x)*cos(y)*sin(ws))

  if (is.null(Rs)) Rs <- 0.5496*Ra

   press = 101.325
  es <- (cal_es(Tmax) + cal_es(Tmin)) / 2
  ea <- SHum * press / (0.378 * SHum + 0.622)

  # 3. Calculate Rn from Rs
  Rns <- (1-0.23)*Rs
  Rso <- (0.75+2*10^(-5)*Z)*Ra
  Rnl <- 4.903*10^(-9)*((Tmax+273.16)^4+(Tmin+273.16)^4)/2*(0.34-0.14*ea^(0.5))*(1.35*Rs/Rso-0.35)
  Rn <- Rns-Rnl
  Rn[Rn < 0] = 0

  # 4. Calculate PET
  delta <- 4098*cal_es(Tm)/(Tm+237.3)^2
  gama <- 0.000665*101.3*((293-0.0065*Z)/293)^5.26
  PET <- (0.408*delta*Rn+gama*900/(Tm+273)*WS*(es-ea))/(delta+gama*(1+0.34*WS))
  return(PET)
}
