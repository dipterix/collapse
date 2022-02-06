#include <cstdlib>
#include <Rcpp.h>
#include "utils.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace Rcpp;
// using namespace RcppParallel;

struct Collapse : public RcppParallel::Worker
{
  const SEXP x_;
  const SEXP dims_;
  const SEXP dim_cumprod;
  const SEXP keep_;
  const SEXP remain;
  const SEXP re;
  const SEXP dim_re;
  const R_xlen_t ndims;
  const SEXP loc_buf;
  const int ncores;

  double* x_ptr;
  R_xlen_t xlen;
  int* dim_ptr;
  int* dim_cumprod_ptr;
  int* keep_ptr;
  int* remain_ptr;
  R_xlen_t re_len;
  double* re_ptr;
  R_xlen_t block_size;
  R_xlen_t nkeeps;
  R_xlen_t collapse_len;

  Collapse(
    const SEXP x,
    const SEXP dims_x,
    const SEXP dim_cumprod,
    const SEXP keep,
    const SEXP remain,
    const SEXP re,
    const SEXP dim_re,
    const R_xlen_t ndims,
    const SEXP loc_buf,
    const int ncores
  ): x_(x), dims_(dims_x), dim_cumprod(dim_cumprod),
  keep_(keep), remain(remain), re(re), dim_re(dim_re),
  ndims(ndims), loc_buf(loc_buf), ncores(ncores){
    this->x_ptr = REAL(this->x_);
    this->xlen = Rf_xlength(this->x_);
    this->dim_ptr = INTEGER(this->dims_);
    this->dim_cumprod_ptr = INTEGER(this->dim_cumprod);
    this->keep_ptr = INTEGER(this->keep_);
    this->remain_ptr = INTEGER(this->remain);
    this->re_len = Rf_xlength(this->re);
    this->re_ptr = REAL(this->re);
    this->block_size = this->re_len / this->ncores;
    if( this->block_size * this->ncores < this->re_len ){
      this->block_size += 1;
    }
    this->nkeeps = Rf_xlength(this->keep_);
    this->collapse_len = this->xlen / this->re_len;
  }

  void per_thread(std::size_t thread) {
    int* dim_re_ptr = INTEGER(dim_re);

    int* dim_cumprod_ptr2 = this->dim_cumprod_ptr;
    int* margin_idx_ptr;
    int* margin_fct_ptr;
    R_xlen_t margin_rem = 0;
    R_xlen_t x_idx = 0;
    double tmp = 0.0;

    int* loc_buf_ptr = INTEGER(loc_buf) + thread * ndims;
    int* loc_buf_ptr2 = loc_buf_ptr;

    // index
    R_xlen_t jj, kk;

    R_xlen_t start = thread * this->block_size;
    R_xlen_t end = (thread + 1) * this->block_size;
    if( end > re_len ){
      end = re_len;
    }
    if(end <= start){ return; }

    for(R_xlen_t ii = start; ii < end; ii++){

      // calculate location indexes on kept margins
      margin_rem = ii;
      // margin_fct = 1;

      margin_idx_ptr = keep_ptr;
      for(jj = 0; jj < nkeeps; jj++){

        // margin_idx = *(keep_ptr + jj);
        // margin_fct = *(dim_re_ptr + jj);
        margin_fct_ptr = dim_re_ptr + jj;
        loc_buf_ptr2 = loc_buf_ptr + *margin_idx_ptr++;
        *loc_buf_ptr2 = margin_rem % *margin_fct_ptr;
        margin_rem = (margin_rem - *loc_buf_ptr2) / *margin_fct_ptr;

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
          loc_buf_ptr2 = loc_buf_ptr + *margin_idx_ptr;
          *loc_buf_ptr2 = margin_rem % *margin_fct_ptr;
          margin_rem = (margin_rem - *loc_buf_ptr2) / *margin_fct_ptr;

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

  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t thread = begin; thread < end; thread++){
      this->per_thread(thread);
    }

  }

};


// [[Rcpp::export]]
SEXP collapser_real_rcppparallel(SEXP x, SEXP keep) {
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
    // R_xlen_t xlen = Rf_xlength(x_);
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


    int ncores = re_len;
    if(ncores > 250000){
      ncores = 250000;
    }

    SEXP loc_buf = PROTECT(Rf_allocVector(INTSXP, ndims * ncores));

    Collapse collapse(x, dims_, dim_cumprod, keep_, remain, re, dim_re, ndims, loc_buf, ncores);
    parallelFor(0, ncores, collapse);

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
