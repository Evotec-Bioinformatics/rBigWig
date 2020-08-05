#' @importFrom assertthat assert_that
#' @useDynLib rBigWig c_fetch_region c_fetch_region_means
#'
NULL



#' Helper function to count the number of links to a file (calls C)
#' @param path A non-empty character string pointing to an existing file
#' @param chromsome
#' @param start
#' @param end
#'
#' @return A data.frame with either 2 or three columns containing either
#'  `Position` and `Score` or `Start`, `End`, `Score` respectively.
#'
#' @export
fetch_region <- function(path, chromosome, start, end) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start > 0)
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

#' Helper function to count the number of links to a file (calls C)
#' @param path A non-empty character string pointing to an existing file
#' @param chromsome
#' @param start
#' @param end
#'
#' @return A data.frame with either 2 or three columns containing either
#'  `Position` and `Score` or `Start`, `End`, `Score` respectively.
#'
#' @export
fetch_region_means <- function(path, chromosome, start, end, bins = 1) {
  assert_that(is.character(path), length(path) == 1)
  assert_that(file.exists(path))
  assert_that(is.character(chromosome), length(chromosome) == 1)
  assert_that(is.numeric(start), start > 0)
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
  d
}


#x <- rBigWig::fetch_region("/data/projects/targetbcd/Alignments/SE2002_1_CSFP200001243-1a_HT7GYDSXX_L1_1_bwa/SE2002_1_CSFP200001243-1a_HT7GYDSXX_L1_1.bw", "21", 1000, 10000)
