#' Compute synthetic estimator. This is a wrapper function for combine_estimators.
#' 
#'
#'
#' @param theta vector or one-row data frame of candidate estimators.
#' @param boot matrix or data frame of bootstrapped versions of candidate estimators.
#' @param theta_0 name of estimator presumed to be unbiased. 
#' @param ... additional arguments to be passed to \code{combine_estimators}.
#' 
#' @return vector of synthetic estimators
#' @export
#'
#' @import dplyr
#' @importFrom purrr simplify

get_synthetic_estimator <- function(theta, boot, theta_0 = NULL, ...) {
  if (!is.data.frame(theta)) {
    theta_ds <- as_tibble(matrix(theta, 1, length(theta)))
    if (!is.null(names(theta))) colnames(theta_ds) <- names(theta)
  } else theta_ds <- theta
  comb <- combine_estimators(ests = theta_ds,
                     boot_ests = as.matrix(boot),
                     name_0 = theta_0,
                     ...)
  coefs <- comb$b_res %>% group_by(theta_0) %>% nest(coef = c(est, b))
  out_df <- comb$ate_res %>%
    transmute(theta_0,
              synthetic_est = ate,
              estimated_mse = var) %>%
    as_tibble %>%
    inner_join(coefs)
  out_df
}


