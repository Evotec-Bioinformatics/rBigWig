#' @importFrom assertthat assert_that
#' @useDynLib rBigWig c_fetch_region c_fetch_region_means c_fetch_region_stdev c_fetch_region_max c_fetch_region_min c_fetch_region_cov c_fetch_region_sum
#'
NULL



#' Fetch the score-values from a BigWig file in the provided region.
#'
#' @param path A non-empty character string pointing to an existing file
#' @param chromsome The name of the chromosome to fetch from
#' @param start An integer value defining the start position (0-based half-open)
#' @param end An integer value defining the end position (0-based half-open). Must be larger than `start`
#'
#' @return A data.frame with either 2 or three columns containing either
#'  `Position` and `Score` or `Start`, `End`, `Score` respectively.
#'
#' @export
fetch_region <- function(path, chromosome, start, end) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start >= 0)
  assert_that(is.numeric(end), end > start)

  d <- .Call(c_fetch_region,
        path[1],
        chromosome[1],
        as.integer(start),
        as.integer(end))
  if (!is.list(d) || length(d) < 2 || length(d) > 3) {
    return(NULL)
  }

  if (length(d) == 2) {
    names(d) <- c("Position", "Score")
  } else {
    names(d) <- c("Start", "End", "Score")
  }

  as.data.frame(d)
}



#' Calculate statistics on `bins` sub-intervals in the provided region
#' 
#' @param path A non-empty character string pointing to an existing file
#' @param chromosome The name of the chromosome to fetch from
#' @param start An integer value defining the start position (0-based half-open)
#' @param end An integer value defining the end position (0-based half-open). Must be larger than `start`.
#' @param bins The number of sub-bins to generate.
#' @param as_data_frame Return a data.frame with precalculated start and ends of the intervals.
#'
#' @return A numeric vector of length `bins` or a `data.frame` with `bins` rows and the columns
#' `Start`, `End`, `Bin`, and `Score`.
#'
#' @name fetch_region_stats
NULL

#' @describeIn fetch_region_stats Calculates the average in the bins
#' @export
fetch_region_means <- function(path, chromosome, start, end, bins = 1, as_data_frame = FALSE) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start >= 0)
  assert_that(is.numeric(end), end > start)

  d <- .Call(c_fetch_region_means,
             path[1],
             chromosome[1],
             as.integer(start),
             as.integer(end),
             as.integer(bins))
  if (!is.vector(d) || !is.numeric(d) || length(d) != bins) {
    return(NULL)
  }


  if (as_data_frame) {
    bin_borders = seq.int(from = start, to = end, length.out = bins + 1)
    data.frame(
      Bin = 1:bins,
      Start = bin_borders[ -length(bin_borders)],
      End = bin_borders[ -1],
      Score = d)
  } else {
    d
  }
}

#' @describeIn fetch_region_stats Calculates the standard-deviation in the bins
#' @export
fetch_region_stdevs <- function(path, chromosome, start, end, bins = 1, as_data_frame = FALSE) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start >= 0)
  assert_that(is.numeric(end), end > start)

  d <- .Call(c_fetch_region_stdev,
             path[1],
             chromosome[1],
             as.integer(start),
             as.integer(end),
             as.integer(bins))
  if (!is.vector(d) || !is.numeric(d) || length(d) != bins) {
    return(NULL)
  }
	if (as_data_frame) {
    bin_borders = seq.int(from = start, to = end, length.out = bins + 1)
    data.frame(
      Bin = 1:bins,
      Start = bin_borders[ -length(bin_borders)],
      End = bin_borders[ -1],
      Score = d)
  } else {
    d
  }
}

#' @describeIn fetch_region_stats Return the maximum coverage value within each bin
#' @export
fetch_region_max <- function(path, chromosome, start, end, bins = 1, as_data_frame = FALSE) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start >= 0)
  assert_that(is.numeric(end), end > start)

  d <- .Call(c_fetch_region_max,
             path[1],
             chromosome[1],
             as.integer(start),
             as.integer(end),
             as.integer(bins))
  if (!is.vector(d) || !is.numeric(d) || length(d) != bins) {
    return(NULL)
  }
	if (as_data_frame) {
    bin_borders = seq.int(from = start, to = end, length.out = bins + 1)
    data.frame(
      Bin = 1:bins,
      Start = bin_borders[ -length(bin_borders)],
      End = bin_borders[ -1],
      Score = d)
  } else {
    d
  }
}



#' @describeIn fetch_region_stats Return the minimum coverage value within each bin
#' @export
fetch_region_min <- function(path, chromosome, start, end, bins = 1, as_data_frame = FALSE) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start >= 0)
  assert_that(is.numeric(end), end > start)

  d <- .Call(c_fetch_region_min,
             path[1],
             chromosome[1],
             as.integer(start),
             as.integer(end),
             as.integer(bins))
  if (!is.vector(d) || !is.numeric(d) || length(d) != bins) {
    return(NULL)
  }
	if (as_data_frame) {
    bin_borders = seq.int(from = start, to = end, length.out = bins + 1)
    data.frame(
      Bin = 1:bins,
      Start = bin_borders[ -length(bin_borders)],
      End = bin_borders[ -1],
      Score = d)
  } else {
    d
  }
}


#' @describeIn fetch_region_stats Return the number of bases covered in each bin
#' @export
fetch_region_covered <- function(path, chromosome, start, end, bins = 1, as_data_frame = FALSE) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start >= 0)
  assert_that(is.numeric(end), end > start)

  d <- .Call(c_fetch_region_cov,
             path[1],
             chromosome[1],
             as.integer(start),
             as.integer(end),
             as.integer(bins))
  if (!is.vector(d) || !is.numeric(d) || length(d) != bins) {
    return(NULL)
  }
	if (as_data_frame) {
    bin_borders = seq.int(from = start, to = end, length.out = bins + 1)
    data.frame(
      Bin = 1:bins,
      Start = bin_borders[ -length(bin_borders)],
      End = bin_borders[ -1],
      Score = d)
  } else {
    d
  }
}


#' @describeIn fetch_region_stats Return the number of bases covered in each bin
#' @export
fetch_region_sum <- function(path, chromosome, start, end, bins = 1, as_data_frame = FALSE) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start >= 0)
  assert_that(is.numeric(end), end > start)

  d <- .Call(c_fetch_region_sum,
             path[1],
             chromosome[1],
             as.integer(start),
             as.integer(end),
             as.integer(bins))
  if (!is.vector(d) || !is.numeric(d) || length(d) != bins) {
    return(NULL)
  }
	if (as_data_frame) {
    bin_borders = seq.int(from = start, to = end, length.out = bins + 1)
    data.frame(
      Bin = 1:bins,
      Start = bin_borders[ -length(bin_borders)],
      End = bin_borders[ -1],
      Score = d)
  } else {
    d
  }
}

