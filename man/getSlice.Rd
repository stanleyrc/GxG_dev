% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apps.R
\name{getSlice}
\alias{getSlice}
\title{getSlice}
\usage{
getSlice(anchors, file, res, gr = NULL, mc.cores = 1)
}
\arguments{
\item{anchors}{result from getBins}

\item{file}{path to mcool}

\item{res}{resolution to use that is in the mcool file}

\item{gr}{optional granges to get from mcool file}

\item{mc.cores}{optional increase the cores when extracting for multiple granges}
}
\value{
list of resolutions in cooler
}
\description{
get resolutions from a mcool file
}
\author{
Stanley Clarke
}
\keyword{getSlice}
