#include "openmp.h"
#include "utils.h"
using namespace Rcpp;



// [[Rcpp::export]]
SEXP collapser_real_omp(SEXP x, SEXP keep) {
  SEXP re = R_NilValue;
  SEXP x_ = R_NilValue;
  if(TYPEOF(x) != REALSXP){
    PROTECT(x_ = Rf_coerceVector(x, REALSXP));
  } else {
    PROTECT(x_ = x);
  }

  R_xlen_t nkeeps = Rf_xlength(keep);
  SEXP keep_ = R_NilValue;
  int* keep_ptr;
  if(TYPEOF(keep) != INTSXP){
    PROTECT(keep_ = Rf_coerceVector(keep, INTSXP));
    keep_ptr = INTEGER(keep_);
    for(R_xlen_t ii = 0; ii < nkeeps; ii++, keep_ptr++){
      *keep_ptr -= 1;
    }
  } else {
    PROTECT(keep_ = Rf_allocVector(INTSXP, nkeeps));
    keep_ptr = INTEGER(keep_);
    for(R_xlen_t ii = 0; ii < nkeeps; ii++){
      *(INTEGER(keep_) + ii) = *(INTEGER(keep) + ii) - 1;
    }
  }
  keep_ptr = INTEGER(keep_);

  // Get dim(x)
  SEXP dims = PROTECT(Rf_getAttrib(x, R_DimSymbol));
  SEXP dims_ = R_NilValue;
  if(TYPEOF(dims) != INTSXP) {
    PROTECT(dims_ = Rf_coerceVector(dims, INTSXP));
  } else {
    PROTECT(dims_ = dims);
  }
  R_xlen_t ndims = Rf_xlength(dims_);
  int* dim_ptr = INTEGER(dims_);


  if( ndims < 2 ){
    re = make_error("x must be an array with at least two margins.");
  } else if(nkeeps < 1 || nkeeps >= ndims){
    re = make_error("`keep` must be a integer vector of positive length, but less than dimension of `x`.");
  } else {
    keep_ptr = INTEGER(keep_);
    for(R_xlen_t ii = 0; ii < nkeeps; ii++, keep_ptr++){
      if(*keep_ptr < 0 || *keep_ptr >= ndims){
        re = make_error("`keep` must be a integer vector of positive length; `keep` cannot exceed maximum dimensions.");
        break;
      }
    }
  }

  if(re == R_NilValue){
    R_xlen_t xlen = Rf_xlength(x_);
    SEXP remain = PROTECT(Rf_allocVector(INTSXP, ndims - nkeeps));
    int* remain_ptr = INTEGER(remain);
    keep_ptr = INTEGER(keep_);
    bool is_in = false;
    for(R_xlen_t ii = 0; ii < ndims; ii++){
      is_in = false;
      for(R_xlen_t jj = 0; jj < nkeeps; jj++){
        if(*(keep_ptr + jj) == ii){
          is_in = true;
          break;
        }
      }
      if(!is_in){
        *remain_ptr++ = ii;
      }
    }


    SEXP dim_cumprod = PROTECT(Rf_allocVector(INTSXP, ndims));

    int* dim_cumprod_ptr = INTEGER(dim_cumprod);
    *dim_cumprod_ptr = 1;
    dim_ptr = INTEGER(dims_);

    for(R_xlen_t ii = 1; ii < ndims; ii++){
      *(dim_cumprod_ptr + ii) = *(dim_cumprod_ptr + (ii-1)) * *(dim_ptr + (ii-1));
    }

    // calculate returned length
    R_xlen_t re_len = 1;
    SEXP dim_re = PROTECT(Rf_allocVector(INTSXP, nkeeps));
    int* dim_re_ptr = INTEGER(dim_re);
    dim_ptr = INTEGER(dims_);
    keep_ptr = INTEGER(keep_);
    for(R_xlen_t ii = 0; ii < nkeeps; ii++, dim_re_ptr++){
      *dim_re_ptr = *(*(keep_ptr + ii) + dim_ptr);
      re_len *= *dim_re_ptr;
    }
    R_xlen_t collapse_len = xlen / re_len;

    PROTECT(re = Rf_allocVector(REALSXP, re_len));
    // double* re_ptr = REAL(re);
    // for(R_xlen_t ii = 0; ii < re_len; ii++){
    //   *re_ptr++ = 0.0;
    // }
    if(nkeeps > 1){
      // R 4.1 seems to allow to set dims even when ndims == 1
      // but just in case
      Rf_setAttrib(re, R_DimSymbol, dim_re);
    }


    int ncores = getThreads();
    if(ncores > re_len){
      ncores = re_len;
    }

    // Allocate buffers
    SEXP loc_buf = PROTECT(Rf_allocVector(INTSXP, ndims * ncores));
    remain_ptr = INTEGER(remain);

#pragma omp parallel num_threads(ncores)
{
  int* keep_ptr = INTEGER(keep_);
  int* dim_re_ptr = INTEGER(dim_re);
  int* dim_ptr = INTEGER(dims_);
  const double* x_ptr = REAL(x_);
  double* re_ptr = REAL(re);
  int* dim_cumprod_ptr2 = dim_cumprod_ptr;
  R_xlen_t margin_rem = 0;
  R_xlen_t margin_rem2 = 0;
  // R_xlen_t margin_fct = 1;
  // R_xlen_t margin_idx = 0;
  int* margin_idx_ptr;
  int* margin_fct_ptr;
  R_xlen_t x_idx = 0;
  double tmp = 0.0;
  int thread = 0;

  // index
  R_xlen_t jj, kk;

#pragma omp for schedule(static, 1) nowait
    for(R_xlen_t ii = 0; ii < re_len; ii++){

      // Get current thread number
      thread = ii % ncores;
      int* loc_buf_ptr = INTEGER(loc_buf) + thread * ndims;
      int* loc_buf_ptr2 = loc_buf_ptr;

      // calculate location indexes on kept margins

      margin_rem = ii;
      // margin_fct = 1;

      margin_idx_ptr = keep_ptr;
      for(jj = 0; jj < nkeeps; jj++){

        // margin_idx = *(keep_ptr + jj);
        // margin_fct = *(dim_re_ptr + jj);
        margin_fct_ptr = dim_re_ptr + jj;
        // loc_buf_ptr2 = loc_buf_ptr + *margin_idx_ptr++;
        // *loc_buf_ptr2 = margin_rem % *margin_fct_ptr;
        // margin_rem = (margin_rem - *loc_buf_ptr2) / *margin_fct_ptr;
        margin_rem2 = margin_rem / (R_xlen_t)(*margin_fct_ptr);
        *(loc_buf_ptr + *margin_idx_ptr++) = margin_rem - margin_rem2 * *margin_fct_ptr;
        margin_rem = margin_rem2;


      }

      tmp = 0.0;


      for(kk = 0; kk < collapse_len; kk++){
        margin_rem = kk;
        // margin_fct = 1;

        for(jj = 0; jj < ndims - nkeeps; jj++){

          // margin_idx = *(remain_ptr + jj);
          // margin_fct = *(dim_ptr + margin_idx);
          margin_idx_ptr = remain_ptr + jj;
          margin_fct_ptr = dim_ptr + *margin_idx_ptr;
          // loc_buf_ptr2 = loc_buf_ptr + *margin_idx_ptr;
          // *loc_buf_ptr2 = margin_rem % *margin_fct_ptr;
          // margin_rem = (margin_rem - *loc_buf_ptr2) / *margin_fct_ptr;

          margin_rem2 = margin_rem / (R_xlen_t)(*margin_fct_ptr);
          *(loc_buf_ptr + *margin_idx_ptr) = margin_rem - margin_rem2 * *margin_fct_ptr;
          margin_rem = margin_rem2;

        }

        // dim_cumprod_ptr = INTEGER(dim_cumprod);
        // loc_buf_ptr = INTEGER(loc_buf);

        loc_buf_ptr2 = loc_buf_ptr;
        dim_cumprod_ptr2 = dim_cumprod_ptr;
        x_idx = 0;
        for(jj = 0; jj < ndims; jj++){
          // print(loc_buf);
          x_idx += *(loc_buf_ptr2++) * *(dim_cumprod_ptr2++);
        }

        tmp += *(x_ptr + x_idx);
        // Rcout << x_idx << "--------\n";
      }

      *(re_ptr + ii) = tmp;

    }
}

    UNPROTECT(5); // remain, re, dim_re_ptr, loc_buf
  }

  UNPROTECT(4); // x_, keep_, dims, dims_
  return re;
//
//   // Generate template output
//   int len = 1;
//   for(int i=0; i<keep.length(); i++){
//     len *= dims[keep[i]-1];
//   }
//   Rcpp::NumericVector re(len);
//
//   // Calculate total dim
//   int total_dim = 1;
//   for(int i=0; i<dims.length(); i++){
//     total_dim *= dims[i];
//   }
//
//   // Calculate remaining dimensions
//   int remsize = dims.size() - keep.size();
//   bool is_in;
//   Rcpp::IntegerVector remain(remsize);
//   for(int64_t j = dims.size(); j > 0; j-- ){
//     is_in = std::find(keep.begin(), keep.end(), j) != keep.end();
//     if(!is_in){
//       remain[--remsize] = j - 1;
//     }
//   }
//
//
//   Collapse collapse(x, dims, keep, remain, total_dim, len, re);
//
//   parallelFor(0, len, collapse);
//
//   return(re);
//   // return(Rcpp::as<Rcpp::NumericVector>(Rcpp::NumericVector::create(len)));
}

/*** R
# x <- matrix(1:16, 4)
# collapse:::baseline_array(x, along_dim = 2, baseline_indexpoints = 1:2)
# dipsaus::baseline_array(x, along_dim = 2, baseline_indexpoints = 1:2)
RcppParallel::setThreadOptions(numThreads = 4)
dat = array(1:16, c(4,4))
dat[1,1] = NA
dat[2,1] = Inf
dat[3,1] = NaN
re = collapser(dat, 1); re
rowSums(dat)
*/
