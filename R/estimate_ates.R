#' Estimate many average treatment effects
#'
#'
#' @param this_data data frame to estimate from.
#' @param fns list of functions to compute ATEs.
#' @param profile logical indicating whether the computation times should be returned rather than the estimates themselves.
#' @param ... additional arguments to be passed to the functions in \code{fns}
#'
#' @return a data frame containing all the ATEs, one per column.
#' @export
#' @importFrom purrr map
#'
#' @examples
#' 
#' 
estimate_ates <- function(this_data, fns, profile = FALSE, ...) {
  # browser()
  ests <- map(fns, function(f) {
    f(this_data, ...)
  }) %>% bind_cols
  if (!profile) {
    return(ests)
  } else {
    times <- sapply(fns, function(f) {
      system.time(f(this_data, ...))[3]
    })
    names(times) <- names(fns)
    return(list(ests = ests, times = times))
  }
}

