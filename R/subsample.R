subsample_1 <- function(data, m, fns, ...) {
  # browser()
  ind <- 1:nrow(data)
  new_ind <- sample(ind, size = m)
  new_data <- data[new_ind,]
  if (all(new_data$d == 0) | all(new_data$d == 1)) return(lapply(fns, function(f) return(NA)) %>% bind_cols)
  new_data <- estimate_scores(new_data, ...)
  lapply(fns, function(f) {
    f(new_data, ...)
  }) %>% bind_cols
}

subsample <- function(data, B = NULL, m = NULL, fns, ...) {
  # browser()
  n <- nrow(data)
  if (is.null(m)) m <- 0.75*n %>% round
  if (is.null(B)) B <- choose(n, m)
  thetahat <- estimate_ates(data, fns, ...) %>% unlist
  ss_fn <- function(x) subsample_1(data, m, fns, ...)
  estimates <- replicate(B, ss_fn(), simplify = FALSE) %>% bind_rows
  # est_cov <- cov(t(t(estimates) - thetahat))
  naive_cov <- cov(estimates, use = 'p')*(m)/(n-m)
  # smart_cov <- (t(estimates) - thetahat) %*% t(t(estimates) - thetahat)/B*(m)/(n-m)
  list(estimates = estimates,
       naive_cov = naive_cov,
       # smart_cov = smart_cov,
       B = B,
       m = m)
}
