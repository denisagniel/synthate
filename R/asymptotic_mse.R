#' Estimate asymptotic MSE
#'
#'
#' @param V_0 scalar variance of theta_0
#' @param V covariance matrix of other candidate estimators
#' @param C vector of covariances of theta_0 with candidates
#' @param deltahat estimate of biases of candidates
#'
#' @return estimated MSE of synthetic estimator
#' @export
#' @importFrom mvnfast rmvn
#'
#' @examples
#' 
#' k <- 5
#' n <- 100
#' theta <- c(0, runif(k-1))
#' Sigma <- 10*diag(k)/n

#' d <- theta[-1]
#' V_0 <- Sigma[1,1]
#' V <- Sigma[-1,-1]
#' C <- Sigma[-1,1]
#' asymptotic_mse(V_0, V, C, d)
#' 
asymptotic_mse <- function(V_0, V, C, deltahat) {
  d <- deltahat
  P <- V_0 - C
  k <- length(C)
  TT <- V_0 + V - (matrix(1,nrow = k) %*% matrix(C, ncol = k)) - 
    (matrix(C,nrow = k) %*% matrix(1, ncol = k))
  J <- mvnfast::rmvn(1000, mu = d, sigma = TT)
  K <- apply(J, 1, function(j) {
    (t(P) %*% solve(TT) %*% j %*% t(j) %*% solve(TT) %*% j) /
      (1 + t(j) %*% solve(TT) %*% j)
  })
  E_K <- mean(K)
  V_K <- var(K)
  E_K2 <- mean(K^2)
  rho <- E_K/(t(P) %*% solve(TT) %*% d)
  lambda <- E_K2/(t(P) %*% solve(TT) %*% (TT + d %*% t(d)) %*% solve(TT) %*% P)
  
  mse_G <- V_0 + (lambda - 1) * t(P) %*% solve(TT) %*% P + (t(P) %*% solve(TT) %*% d)^2 * (2*(1-rho) - 1 + lambda)
  mse_G
}
