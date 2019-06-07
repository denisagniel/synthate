#' Generate data for compliance using the Ding and Lu model. 
#'
#' @param n sample size
#' @param theta paramater controlling the violation of principal ignorability
#' @return data.frame including \code{id}, covariate \code{x}, randomization indicator \code{z}, complier indicator \code{c}, treatment receipt indicator \code{s}, and outcome \code{y}.
#' @export
#'
#' @examples
#' generate_data_dl(n = 500, theta = 1)

generate_data_dl <- function(n = 500, 
                             theta) {
  # browser()
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  x5 <- rbinom(n, size = 1, prob = 0.5)
  
  pi_c <- plogis(0.5*x1 + 0.5*x2 + x3 + x4 + theta*x5)
  
  z <- rbinom(n, prob = 0.5, size = 1) ## randomization indicator
  c <- rbinom(n, prob = pi_c, size = 1) ## compliance indicator
  epsilon <- rnorm(n)
  
  
  y_1 <- x1 + x2 + x3 + x4 + x5 + 2*c + 1 + epsilon
  y_0 <- x1 + x2 + x3 + x4 + x5 + 2 + epsilon
  y <- if_else(z == 1 & c == 1, y_1, y_0)
  
  sim_data <- tibble(id = 1:n,
                         x = cbind(x1, x2, x3, x4),
                         z = z,
                         c = c,
                         s = 1*(z == 1 & c == 1), ## receipt of treatment
                         y = y)
  sim_data
}
