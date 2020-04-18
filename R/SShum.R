#' calculate saturated specific humidity for MEP model input
#' @param  Ts surface temperature(unit: deg C)
#' @return saturated specific humidity(unit:kg/kg)
#' @examples
#' SShum(Ts=20)
#' @export
SShum <- function(Ts){
  t_temp<-C2K(Ts)
  Es<-SVP.ClaCla(t_temp)  #calculate saturated vapor pressure Es
  Sqs<-SH(Es*100)        #calculate saturated specific humidity
  return(Sqs)
}
