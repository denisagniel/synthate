#' Generate data for compliance under pseudo-randomization, where pseudo-randomization holds within strata.
#'
#' @param n sample size
#' @param p dimension of covariates
#' @param delta complier average causal effect
#' @param compliance_p compliance proportion
#'
#' @return data.frame including \code{id}, treatment assignment indicator \code{z}, complier indicator \code{c}, treatment receipt indicator \code{s}, outcome \code{y}, and covariates \code{x1, ..., xp}.
#' @importFrom withr with_seed
#' @importFrom tidyr spread
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' generate_data_psrg(n = 100, p = 8)

generate_data_psrg <- function(n, 
                              p,
                              delta = 0.5, 
                              compliance_p = 0.75
) {
  # browser()
  x <- matrix(rnorm(n*p), n, p) ## covariates
  # w <- rnorm(n) ## additional covariate that only determines compliance
  theta <- withr::with_seed(0, rnorm(p))
  eta <- withr::with_seed(1, rnorm(p))
  gamma_1 <- withr::with_seed(2, rnorm(p))
  gamma_0 <- withr::with_seed(3, rnorm(p))
  
  p_z <- plogis(x %*% theta)
  g_z <- Hmisc::cut2(p_z, g = 5)
  pz_ds <- tibble::tibble(p_z, g_z)
  pz_ds <- dplyr::group_by(pz_ds, g_z)
  pz_ds <- mutate(pz_ds, pi_z = mean(p_z))
  pi_z <- dplyr::pull(pz_ds, pi_z)
  pi_c <- plogis(x %*% eta + qlogis(compliance_p)) ## probability of compliance
  z <- rbinom(n, prob = pi_z, size = 1) ## assignment indicator
  c <- rbinom(n, prob = pi_c, size = 1) ## compliance indicator
  epsilon <- rnorm(n)
  
  
  y_1c <- -1 + x %*% gamma_1 + delta + epsilon
  y_0c <- -1 + x %*% gamma_1 + epsilon
  y_1n <- x %*% gamma_0 + epsilon
  y_0n <- x %*% gamma_0 + epsilon
  
  sim_data <- tibble::tibble(id = rep(1:n, p),
                             x = c(x),
                             x_n = rep(paste('x', 1:p, sep = ''), each = n),
                             z = rep(z, p),
                             c = rep(c, p),
                             s = 1*(z == 1 & c == 1), ## receipt of treatment
                             y1c = rep(y_1c,p),
                             y0c = rep(y_0c,p),
                             y1n = rep(y_1n,p),
                             y0n = rep(y_0n,p),
                             y = dplyr::case_when(
                               s == 1 & c == 1 ~ y1c,
                               s == 0 & c == 1 ~ y0c,
                               s == 1 & c == 0 ~ y1n,
                               s == 0 & c == 0 ~ y0n)
  )
  sim_data <- tidyr::spread(sim_data, x_n, x)
  sim_data
}
