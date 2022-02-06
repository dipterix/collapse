
stopifnot2 <- function(..., msg = 'Condition not satisfied'){
  if(!all(c(...))){
    stop(msg)
  }
}

.missing_arg <- alist(x = )
