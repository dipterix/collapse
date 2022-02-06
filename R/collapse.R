# baseline_array <- function(x, along_dim, baseline_indexpoints,
#                            unit_dims = seq_along(dim(x))[-along_dim],
#                            method = c("percentage", "sqrt_percentage",
#                                       "decibel", "zscore", "sqrt_zscore")) {
#   # baselineArray(x, bl, dim(x), dim(bl), 2L, c(0L, 1L, 3L), c(2L), method)
#   stopifnot2(is.array(x), msg = paste0(sQuote('x'), ' must be an array/matrix'))
#   dims <- dim(x)
#   along_dim <- as.integer(along_dim)
#   stopifnot2(along_dim >=1 && along_dim <= length(dims), msg = paste0(
#     sQuote('along_dim'), ' is invalid, it must be an integer from 1 to ',
#     length(dims)
#   ))
#
#   # method
#   blmethods <- c("percentage", "sqrt_percentage", "decibel", "zscore", "sqrt_zscore")
#   method <- match.arg(method)
#   method <- which(blmethods %in% method)
#   stopifnot2(length(method) == 1, msg = paste0(
#     sQuote('method'), ' must be one of the following options: ',
#     paste(blmethods, collapse = ', ')
#   ))
#
#   # unit_dims is baseline unit,
#   unit_dims <- as.integer(unit_dims)
#   stopifnot2(!any(is.na(unit_dims)), msg = paste0(sQuote('unit_dims'), ' contains NAs'))
#   stopifnot2(all(unit_dims %in% seq_along(dims)), msg = paste0(sQuote('unit_dims'), ' has invalid dimensions'))
#   stopifnot2(!along_dim %in% unit_dims, msg = paste0(
#     sQuote('along_dim'), ' cannot be inside of ', sQuote('unit_dims')))
#   unit_dims <- sort(unit_dims)
#
#   ntimepoints <- dims[along_dim]
#   # calculate baseline window
#   baseline_indexpoints <- as.integer(baseline_indexpoints)
#   if(any(is.na(baseline_indexpoints))){
#     baseline_indexpoints <- baseline_indexpoints[!is.na(baseline_indexpoints)]
#   }
#   if(!is.unsorted(baseline_indexpoints)){
#     baseline_indexpoints <- sort(baseline_indexpoints)
#   }
#   baseline_indexpoints <- baseline_indexpoints[baseline_indexpoints > 0 & baseline_indexpoints <= ntimepoints]
#   stopifnot2(length(baseline_indexpoints) > 0, msg = paste0('Baseline window is invalid: cannot find any valid time points. \nPlease makesure ', sQuote('baseline_indexpoints'), ' contains at least one integer from 1 to ', ntimepoints))
#
#   call_args <- lapply(seq_along(dims), function(ii){
#     if(ii == along_dim){
#       return(quote(baseline_indexpoints))
#     }
#     .missing_arg[[1]]
#   })
#   call_args <- c(list(quote(x)), call_args, list(drop = FALSE))
#   bl <- do.call('[', call_args)
#   bldims <- dim(bl)
#
#   rest <- seq_along(dims)[-unit_dims]
#
#   .Call(`_collapse_baselineArray`, x, bl, dims, bldims, along_dim - 1L, unit_dims - 1L, rest - 1L, method - 1L)
# }
