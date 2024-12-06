% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strip_transect.R
\name{strip_transect}
\alias{strip_transect}
\title{Create a strip map with raster data using transects. EXPERIMENTAL}
\usage{
strip_transect(line, r, shape, n_transects, width, radius, ...)
}
\arguments{
\item{line}{An sf linestring}

\item{r}{sf object with the vectors to be strip'd.}

\item{shape}{a sf linestring object showing the target shape.}

\item{n_transects}{The number of transects to create from the line}

\item{width}{Width of the transect}

\item{radius}{radius for inverse distance weighting-based interpolation to
fill gaps between transects.}

\item{...}{Further arguments passed to terra::interpIDW()}
}
\value{
a list with two elements: the rotated linestring and vectors.
}
\description{
The function creates a strip of the specified shape with raster data 
using transects derived from the linestring. Raster data is then sampled
from the raster with the linestring. The values are then translated to
transects computed from the target shape, and gaps are filled with 
using IDW interpolation. EXPERIMENTAL AND SUBJECT TO CHANGE!
}
