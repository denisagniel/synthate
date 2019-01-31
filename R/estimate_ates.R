estimate_ates <- function(this_data, fns, profile = FALSE, ...) {
  # browser()
  ests <- lapply(fns, function(f) {
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

