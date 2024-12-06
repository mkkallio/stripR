% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/straighten_line.R
\name{straighten_line}
\alias{straighten_line}
\title{Straighten a linestring to create a strip}
\usage{
straighten_line(sf, tolerance, angle = 0)
}
\arguments{
\item{sf}{An sf linestring}

\item{tolerance}{tolerance input for the Douglas-Peucker algorithm.}

\item{angle}{The direction of of the straightened line. Default = 0 (north).}
}
\value{
a list with 1. the straightened line, and 2. the original line with
segments defined by the simplification.
}
\description{
The function creates a straight strip from the linestring provided. Douglas-
Peucker simplification algorithm is used to find nodes which are used to 
split the linestring in to segments. The segments are then rotated to
the angle provided, and stitched together. The result is a straight line
with details of the original linestring. The width of the strip (and amount
of detail in the linestring) is controlled by the tolerance argument.
}
