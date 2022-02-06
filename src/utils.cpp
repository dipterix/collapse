#include "utils.h"

SEXP make_error(const char* message){
  SEXP error;
  PROTECT(error = Rf_mkString(message));
  Rf_classgets(error, Rf_mkString("collapse_error"));
  UNPROTECT(1);
  return error;
}

