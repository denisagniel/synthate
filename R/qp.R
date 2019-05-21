#' Solve quadratic program
#'
#'
#' @param qq objective function
#' @param n_ests number of candidate estimators
#'
#' @return list of objects from quadprog::solve.QP
#' @importFrom matrixcalc is.positive.definite
#' @importFrom quadprog solve.QP
#'
#' @examples
#' rho <- 0.5
#' d <- 3
#' cov_mat <- rho*matrix(1, d, d) + (1-rho)*diag(d)
#' bias <- runif(d)
#' mse <- cov_mat + bias %*% t(bias)
#' qp(mse, d)
qp <- function(qq, n_ests) {
  if (!matrixcalc::is.positive.definite(round(qq, 9))) {
    qq_tst <- try(Matrix::nearPD(round(qq, 9))$mat)
    if (class(qq_tst) == 'try-error') browser()
    qq <- qq_tst
    # qq <- Matrix::nearPD(C)$mat %>% as.matrix
  }
  Amat <- cbind(rep(1, n_ests), diag(n_ests))
  quadprog::solve.QP(qq, dvec = rep(0, n_ests),
                     Amat = Amat, bvec = c(1, rep(0, n_ests)),
                     meq = 1)
}
