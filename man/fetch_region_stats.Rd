% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rBigWig.R
\name{fetch_region_stats}
\alias{fetch_region_stats}
\alias{fetch_region_means}
\alias{fetch_region_stdevs}
\alias{fetch_region_max}
\alias{fetch_region_min}
\alias{fetch_region_covered}
\alias{fetch_region_sum}
\title{Calculate statistics on `bins` sub-intervals in the provided region}
\usage{
fetch_region_means(
  path,
  chromosome,
  start,
  end,
  bins = 1,
  as_data_frame = FALSE
)

fetch_region_stdevs(
  path,
  chromosome,
  start,
  end,
  bins = 1,
  as_data_frame = FALSE
)

fetch_region_max(path, chromosome, start, end, bins = 1, as_data_frame = FALSE)

fetch_region_min(path, chromosome, start, end, bins = 1, as_data_frame = FALSE)

fetch_region_covered(
  path,
  chromosome,
  start,
  end,
  bins = 1,
  as_data_frame = FALSE
)

fetch_region_sum(path, chromosome, start, end, bins = 1, as_data_frame = FALSE)
}
\arguments{
\item{path}{A non-empty character string pointing to an existing file}

\item{chromosome}{The name of the chromosome to fetch from}

\item{start}{An integer value defining the start position (0-based half-open)}

\item{end}{An integer value defining the end position (0-based half-open). Must be larger than `start`.}

\item{bins}{The number of sub-bins to generate.}

\item{as_data_frame}{Return a data.frame with precalculated start and ends of the intervals.}
}
\value{
A numeric vector of length `bins` or a `data.frame` with `bins` rows and the columns
`Start`, `End`, `Bin`, and `Score`.
}
\description{
Calculate statistics on `bins` sub-intervals in the provided region
}
\section{Functions}{
\itemize{
\item \code{fetch_region_means}: Calculates the average in the bins

\item \code{fetch_region_stdevs}: Calculates the standard-deviation in the bins

\item \code{fetch_region_max}: Return the maximum coverage value within each bin

\item \code{fetch_region_min}: Return the minimum coverage value within each bin

\item \code{fetch_region_covered}: Return the number of bases covered in each bin

\item \code{fetch_region_sum}: Return the number of bases covered in each bin
}}

