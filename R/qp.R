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
