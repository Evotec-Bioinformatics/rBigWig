#include <sys/types.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <stdio.h>
#include <unistd.h>
#include <R.h>
#include <Rinternals.h>

#include "libBigWig/bigWig.h"

SEXP c_fetch_region(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend) {
  int nprotect = 0;
  PROTECT(Rfilename); nprotect++;
  char* pathname = CHAR(asChar(Rfilename));

  PROTECT(Rchromosome); nprotect++;
  char* chromosome = CHAR(asChar(Rchromosome));

  int start = asInteger(Rstart);
  int end = asInteger(Rend);


  bigWigFile_t *fp = NULL;
  bwOverlappingIntervals_t *intervals = NULL;
  double *stats = NULL;

  //Initialize enough space to hold 128KiB (1<<17) of data at a time
  if(bwInit(1<<17) != 0) {
    fprintf(stderr, "Received an error in bwInit\n");
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Open the local/remote file
  fp = bwOpen(pathname, NULL, "r");
  if(!fp) {
    fprintf(stderr, "An error occured while opening %s\n", pathname);
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Get values in a range (0-based, half open) without NAs
  intervals = bwGetValues(fp, chromosome, start, end, 0);
  int n_intervals = intervals->l;

  SEXP starts = PROTECT(allocVector(INTSXP, n_intervals)); nprotect++;
  if (isNull(starts)) {
    fprintf(stderr, "Can not allocate %i-int vector for starts\n", n_intervals);
    UNPROTECT(nprotect);
    return R_NilValue;
  }
  int* ptx_starts = INTEGER(starts);

  SEXP values = PROTECT(allocVector(REALSXP, n_intervals)); nprotect++;
  if (isNull(values)) {
    fprintf(stderr, "Can not allocate %i-double vector for values\n", n_intervals);
    UNPROTECT(nprotect);
    return R_NilValue;
  }
  double* ptx_values = REAL(values);

  bool has_ends = intervals->end > 0;
  int* ptx_ends = NULL;
  SEXP ends = R_NilValue;
  if (has_ends) {
    ends = PROTECT(allocVector(INTSXP, n_intervals)); nprotect++;
    if (isNull(ends)) {
      fprintf(stderr, "Can not allocate %i-int vector for ends\n", n_intervals);
      UNPROTECT(nprotect);
      return R_NilValue;
    }
    ptx_ends = INTEGER(ends);
  }

  for (int i = 0; i < n_intervals; i++) {
    ptx_starts[i] = intervals->start[i];
    ptx_values[i] = intervals->value[i];

    if (has_ends) {
      ptx_ends[i] = intervals->end[i];
    }
  }

  bwDestroyOverlappingIntervals(intervals); //Free allocated memory

  bwClose(fp);
  bwCleanup();

  SEXP res = R_NilValue;
  if (has_ends) {
    res = PROTECT(allocVector(VECSXP, 3)); nprotect++;
    SET_VECTOR_ELT(res, 0, starts);
    SET_VECTOR_ELT(res, 1, ends);
    SET_VECTOR_ELT(res, 2, values);
  } else {
    res = PROTECT(allocVector(VECSXP, 2)); nprotect++;
    SET_VECTOR_ELT(res, 0, starts);
    SET_VECTOR_ELT(res, 1, values);
  }

  UNPROTECT(nprotect);
  return res;
}


SEXP c_fetch_region_means(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins) {
  int nprotect = 0;
  PROTECT(Rfilename); nprotect++;
  char* pathname = CHAR(asChar(Rfilename));

  PROTECT(Rchromosome); nprotect++;
  char* chromosome = CHAR(asChar(Rchromosome));

  int start = asInteger(Rstart);
  int end = asInteger(Rend);
  int bins = asInteger(Rbins);


  bigWigFile_t *fp = NULL;
  bwOverlappingIntervals_t *intervals = NULL;
  double *stats = NULL;

  //Initialize enough space to hold 128KiB (1<<17) of data at a time
  if(bwInit(1<<17) != 0) {
    fprintf(stderr, "Received an error in bwInit\n");
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Open the local/remote file
  fp = bwOpen(pathname, NULL, "r");
  if(!fp) {
    fprintf(stderr, "An error occured while opening %s\n", pathname);
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Get values in a range (0-based, half open) without NAs
  stats = bwStats(fp, chromosome, start, end, bins, mean);
  if (stats == NULL) {
    fprintf(stderr, "Failed to generate stats in bwStats\n");
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  SEXP values = PROTECT(allocVector(REALSXP, bins)); nprotect++;
  if (isNull(values)) {
    fprintf(stderr, "Can not allocate %i-double vector to store the stats\n", bins);
    UNPROTECT(nprotect);
    return R_NilValue;
  }
  double* ptx_values = REAL(values);

  void *ret = memcpy(ptx_values, stats, sizeof(double) * bins);
  if (ret == NULL) {
    fprintf(stderr, "Can not copy %i-double vector to output R-struct\n", bins);
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  bwClose(fp);
  bwCleanup();

  UNPROTECT(nprotect);
  return values;
}
