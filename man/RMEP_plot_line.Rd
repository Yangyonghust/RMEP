% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotfunction.R
\name{RMEP_plot_line}
\alias{RMEP_plot_line}
\title{Line plot between MEP estimations and observations}
\usage{
RMEP_plot_line(timeseries, labx, laby)
}
\arguments{
\item{timeseries}{A dataframe contains three columns which named time series, MEP estimate and observations, respectively.}

\item{labx}{Label for X-axis}

\item{laby}{Label for Y-axis}
}
\description{
Generating a line plot from a dataframe contains three columns: time series, MEP estimate and observations.
}
\examples{
RMEP_plot_line(timeseries)
}
