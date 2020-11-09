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
  const char* pathname = CHAR(asChar(Rfilename));

  PROTECT(Rchromosome); nprotect++;
  const char* chromosome = CHAR(asChar(Rchromosome));

  int start = asInteger(Rstart);
  int end = asInteger(Rend);


  bigWigFile_t *fp = NULL;
  bwOverlappingIntervals_t *intervals = NULL;
  double *stats = NULL;

  //Initialize enough space to hold 128KiB (1<<17) of data at a time
  if(bwInit(1<<17) != 0) {
    //fprintf(stderr, "Received an error in bwInit\n");
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Open the local/remote file
  fp = bwOpen((char *)pathname, NULL, "r");
  if(!fp) {
    //fprintf(stderr, "An error occured while opening %s\n", pathname);
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Get values in a range (0-based, half open) without NAs
  intervals = bwGetValues(fp, (char *)chromosome, start, end, 0);
  int n_intervals = intervals->l;

  SEXP starts = PROTECT(allocVector(INTSXP, n_intervals)); nprotect++;
  if (isNull(starts)) {
    //fprintf(stderr, "Can not allocate %i-int vector for starts\n", n_intervals);
    UNPROTECT(nprotect);
    return R_NilValue;
  }
  int* ptx_starts = INTEGER(starts);

  SEXP values = PROTECT(allocVector(REALSXP, n_intervals)); nprotect++;
  if (isNull(values)) {
    //fprintf(stderr, "Can not allocate %i-double vector for values\n", n_intervals);
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
      //fprintf(stderr, "Can not allocate %i-int vector for ends\n", n_intervals);
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


SEXP c_fetch_region_stats(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins, SEXP Rwant_data_frame, enum bwStatsType statsMode) {
  int nprotect = 0;
  PROTECT(Rfilename); nprotect++;
  const char* pathname = CHAR(asChar(Rfilename));

  PROTECT(Rchromosome); nprotect++;
  const char* chromosome = CHAR(asChar(Rchromosome));

  int start = asInteger(Rstart);
  int end = asInteger(Rend);
  int bins = asInteger(Rbins);
  bool want_data_frame = asInteger(Rwant_data_frame) != 0;

  bigWigFile_t *fp = NULL;
  bwOverlappingIntervals_t *intervals = NULL;
  double *stats = NULL;

  //Initialize enough space to hold 128KiB (1<<17) of data at a time
  if(bwInit(1<<17) != 0) {
    //fprintf(stderr, "Received an error in bwInit\n");
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Open the local/remote file
  fp = bwOpen((char *)pathname, NULL, "r");
  if(!fp) {
    //fprintf(stderr, "An error occured while opening %s\n", pathname);
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  //Get values in a range (0-based, half open) without NAs
  stats = bwStats(fp, (char *)chromosome, start, end, bins, statsMode);
  if (stats == NULL) {
    //fprintf(stderr, "Failed to generate stats in bwStats\n");
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  SEXP values = PROTECT(allocVector(REALSXP, bins)); nprotect++;
  if (isNull(values)) {
    //fprintf(stderr, "Can not allocate %i-double vector to store the stats\n", bins);
    UNPROTECT(nprotect);
    return R_NilValue;
  }
  double* ptx_values = REAL(values);

  void *ret = memcpy(ptx_values, stats, sizeof(double) * bins);
  if (ret == NULL) {
    //fprintf(stderr, "Can not copy %i-double vector to output R-struct\n", bins);
    UNPROTECT(nprotect);
    return R_NilValue;
  }

  bwClose(fp);
  bwCleanup();

	// Check if a data frame must be returned
	if (!want_data_frame) {
  	UNPROTECT(nprotect);
		return values;
	}

	// Calculate the bin-size
	int bin_size = (end - start) / bins;

	// Calculate the bin-sizes/positions
  SEXP df_start = PROTECT(allocVector(INTSXP, bins)); nprotect++;
  SEXP df_end   = PROTECT(allocVector(INTSXP, bins)); nprotect++;
	SEXP rownames = PROTECT(allocVector(INTSXP, bins)); nprotect++; // create a vector for the rownames; needed because rownames length defines the data.frame row number
	int pos = start;
	for (int i = 0, pos=start; i < bins; i++) {
		INTEGER(df_start)[i] = pos;
    pos = start + ((double)(end-start)*(i+1)) / bins;
		INTEGER(df_end)[i] = pos;
		INTEGER(rownames)[i] = i + 1;
	}

	// Generate a list and add the three elements
	SEXP df = PROTECT(Rf_allocVector(VECSXP, 3)); nprotect++; // a list with three elements: start, end, value
	SET_VECTOR_ELT(df, 0, df_start);
	SET_VECTOR_ELT(df, 1, df_end);
	SET_VECTOR_ELT(df, 2, values);

	// Set class of list to data frame
	SEXP cls = PROTECT(allocVector(STRSXP, 1)); nprotect++; // class attribute
  SET_STRING_ELT(cls, 0, mkChar("data.frame"));
  classgets(df, cls);

	// Specify the column names
	SEXP colnames = PROTECT(allocVector(STRSXP, 3)); nprotect++; // names attribute (column names)
  SET_STRING_ELT(colnames, 0, mkChar("Start"));
  SET_STRING_ELT(colnames, 1, mkChar("End"));
  SET_STRING_ELT(colnames, 2, mkChar("Score"));

	setAttrib(df, R_NamesSymbol, colnames);
	setAttrib(df, R_RowNamesSymbol, rownames);

  UNPROTECT(nprotect);
  return df;
}


SEXP c_fetch_region_means(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins, SEXP Rwant_data_frame) {
  return c_fetch_region_stats(Rfilename, Rchromosome, Rstart, Rend, Rbins, Rwant_data_frame, mean);
}

SEXP c_fetch_region_stdev(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins, SEXP Rwant_data_frame) {
  return c_fetch_region_stats(Rfilename, Rchromosome, Rstart, Rend, Rbins, Rwant_data_frame, stdev);
}

SEXP c_fetch_region_max(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins, SEXP Rwant_data_frame) {
  return c_fetch_region_stats(Rfilename, Rchromosome, Rstart, Rend, Rbins, Rwant_data_frame, max);
}

SEXP c_fetch_region_min(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins, SEXP Rwant_data_frame) {
  return c_fetch_region_stats(Rfilename, Rchromosome, Rstart, Rend, Rbins, Rwant_data_frame, min);
}

SEXP c_fetch_region_cov(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins, SEXP Rwant_data_frame) {
  return c_fetch_region_stats(Rfilename, Rchromosome, Rstart, Rend, Rbins, Rwant_data_frame, cov);
}

SEXP c_fetch_region_sum(SEXP Rfilename, SEXP Rchromosome, SEXP Rstart, SEXP Rend, SEXP Rbins, SEXP Rwant_data_frame) {
  return c_fetch_region_stats(Rfilename, Rchromosome, Rstart, Rend, Rbins, Rwant_data_frame, sum);
}

