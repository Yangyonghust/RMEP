#' Line plot between MEP estimations and observations
#'
#' Generating a line plot from a dataframe contains three columns: time series, MEP estimate and observations.
#' @param timeseries A dataframe contains three columns which named time series, MEP estimate and observations, respectively.
#' @param labx Label for X-axis
#' @param laby Label for Y-axis
#' @importFrom ggplot2 ggplot geom_line aes theme_bw labs theme element_blank scale_colour_manual
#' @examples
#' RMEP_plot_line(timeseries)
#' @export
RMEP_plot_line <- function(timeseries,labx,laby){  ### timeseries contains three variables: time, MEP, obs
  labx="time";laby="Flux";
  ggplot(data = timeseries)+geom_line(aes(x=time,y=MEP,color="MEP"),lwd=0.6)+geom_line(aes(x=time,y=obs,color="Obs"),lwd=0.6)+
    theme_bw()+labs(x=labx,y=laby,fill="")+
    theme(legend.title=element_blank(),legend.position = c(0.9, 0.9))+scale_colour_manual(values = c("#3f60aa","red"))
}

#' Scatter plot between MEP estimations and observations
#'
#' Generating a scatter plot with MEP estimations and observations. The labels contains: label plot with the regression equation (a fitted polynomial model) and label with coefficient of determination (R^2) for fitted model.
#' @param MEP A vector of MEP estimations
#' @param obs A vector of observations
#' @param labx Label for X-axis
#' @param laby Label for Y-axis
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom ggplot2 aes  geom_abline labs theme_bw geom_point geom_smooth
#' @examples
#' RMEP_plot_sca(MEP,obs)
#' @export
RMEP_plot_sca <- function(obs,MEP,labx,laby){
  #library(ggpmisc)
  labx="Obs";laby="MEP";
  ggplot(mapping = aes(x=obs,y=MEP))+geom_point(col="blue")+
    geom_smooth(method="lm",formula = y~x,col="red",lty=1)+
    geom_abline(intercept=0, slope=1,col="black",lty=5,lwd=1)+labs(x="MEP",y="observation")+
    labs(x=labx,y=laby,fill="")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = '~~~~~~')),
                 formula=y~x,parse=TRUE,size=5,label.x=0.05,label.y=0.9)+theme_bw()
}

#' Spatial distribution plot
#'
#' Generating a spatial distribution map of evapotranspiration.The plot includes image and legend, generated based on the fields package.
#' @param lon The varaible of longitude
#' @param lat The variable of latitude
#' @param MEP The variable of MEP estimated evapotranspiration
#' @param main The title of figure
#' @param legends The legend of figure
#' @importFrom fields image.plot
#' @examples
#' RMEP_plot_spa(lon,lat,ETannual)
#' @export
RMEP_plot_spa <- function(lon,lat,MEP,main="Latent heat flux by MEP",legends="W m-2"){
  image.plot(lon,lat,MEP,main=main,legend.lab=legends)
}



