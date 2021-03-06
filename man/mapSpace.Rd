% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapSpace.R
\name{mapSpace}
\alias{mapSpace}
\title{Construct species distribution maps models}
\usage{
mapSpace(
  modelSpace,
  dims,
  type = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode", "space"),
  sPoly = NULL
)
}
\arguments{
\item{modelSpace}{An object of class \code{ppSpace} or of class \code{uniSpace}.}

\item{dims}{A vector of length 2 defining the number of pixels to use as rows and columns to define the map.}

\item{type}{Either "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode" or "space". Defines the map to be drawn. If "space", maps the mean of the spatial field.}

\item{sPoly}{A spatial polygon to isolate the region of interest. If none is given, a map is drawn for the entire region covered by the mesh.}
}
\description{
Constructs mean, standard deviation, and quantile (0.025, 0.5 and 0.975) maps for models calculated using \code{\link{uniSpace}} and \code{\link{ppSpace}}
}
\keyword{hplot}
