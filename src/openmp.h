#ifndef COLLAPSE_OPENMP_H
#define COLLAPSE_OPENMP_H

#ifdef _OPENMP
#include <omp.h>
#include <pthread.h>
#define COLLAPSE_HAS_OPENMP true
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define COLLAPSE_HAS_OPENMP false
#endif

#include "common.h"
#include <algorithm>


// [[Rcpp::export]]
int getThreads(bool max = false);

// [[Rcpp::export]]
int setThreads(int n, int reset_after_fork = -1);

// [[Rcpp::export]]
bool hasOpenMP();

// [[Rcpp::init]]
int detectForked(DllInfo *dll);

#endif  // COLLAPSE_OPENMP_H
