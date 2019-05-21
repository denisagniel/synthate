#' Generate data for compliance using the Stuart and Jo model. 
#'
#' @param gamma_c effect of randomization among the compliers
#' @param compliance_p proportion of compliance
#' @param compliance_effect effect of covariate on compliance
#' @param alpha_c intercept term for compliers
#' @param alpha_n intercept term for never-takers
#' @param lambda_c covariate effect among compliers
#' @param lambda_n covariate effect among never-takers
#' @param gamma_n effect of randomization among never-takers
#' @param sigma_c error variance among compliers
#' @param sigma_n error variance among never-takers
#'
#' @return data.frame including \code{id}, covariate \code{x}, randomization indicator \code{z}, complier indicator \code{c}, treatment receipt indicator \code{s}, and outcome \code{y}.
#' @export
#'
#' @examples
#' generate_data_sj(n = 100)

generate_data_sj <- function(n, 
                             gamma_c = 0.5, 
                             compliance_p = 0.5, 
                             compliance_effect = 0.5, 
                             alpha_c = 1, 
                             alpha_n = 0.5, 
                             lambda_c = 1, 
                             lambda_n = 1, 
                             gamma_n = 0, 
                             sigma_c = 1, 
                             sigma_n = 1) {
  # browser()
  x <- rnorm(n) ## covariate
  z <- rbinom(n, prob = 0.5, size = 1) ## randomization indicator
  if (compliance_effect == 0) {
    pi_c <- plogis(qlogis(compliance_p)) ## probability of compliance
  } else {
    pi_c <- plogis(log(compliance_effect)*x + qlogis(compliance_p)) ## probability of compliance
  }
  
  c <- rbinom(n, prob = pi_c, size = 1) ## compliance indicator
  epsilon <- rnorm(n)
  
  
  y <- alpha_n + (alpha_c - alpha_n)*c + gamma_n*z + (gamma_c - gamma_n)*c*z +
    lambda_n*x + (lambda_c - lambda_n)*c*x + sigma_n*epsilon + (sigma_c - sigma_n)*c*epsilon
  
  sim_data <- data.frame(id = 1:n,
                         x = x,
                         z = z,
                         c = c,
                         s = 1*(z == 1 & c == 1), ## receipt of treatment
                         y = y)
  sim_data
}
