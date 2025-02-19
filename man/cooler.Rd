% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apps.R
\name{cooler}
\alias{cooler}
\title{cooler}
\usage{
cooler(file, res = NULL, gr = NULL, info = FALSE, mc.cores = 1)
}
\arguments{
\item{file}{path to .cool or m.cool file}

\item{res}{optional resolution to pull}

\item{gr}{optional GRanges of intervals to query}

\item{info}{logical flag (FALSE) specifying whether to return a vector of file info instead of pulling down data, this will include available resolutions (for .mcool format) and fields}
}
\value{
gMatrix or list
}
\description{
Instantiates gMatrix from .cool and/or .mcool file, +/- at specific locations. 

To find viable resolutions (e.g. for mcool files) use info = TRUE and you will get a list
of info, including seqlengths.
}
\author{
Marcin Imielinski, Stanley Clarke
}
\keyword{cooler}
