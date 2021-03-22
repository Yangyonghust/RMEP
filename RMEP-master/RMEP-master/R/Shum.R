#' calculate specific humidity for MEP model input
#' @param TA Air temperature(unit: deg C)
#' @param RH Relative humidity(unit: percent)
#' @param PA Atmospheric pressure(unit: kPa)
#' @return specific humidity(unit:kg/kg)
#' @examples
#' Shum(TA=20,PA=101,RH=50)
#' @export
Shum <- function(TA,PA,RH){
  Tk <- C2K(TA)
  Es <- SVP(Tk)
  E <- WVP2(RH,Es)
  qs <- SH(E,PA*1000)   #specific humidity in kg/kg
  return(qs)
}

