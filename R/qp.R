qp <- function(qq, n_ests) {
  Amat <- cbind(rep(1, n_ests), diag(n_ests))
  quadprog::solve.QP(qq, dvec = rep(0, n_ests),
                     Amat = Amat, bvec = c(1, rep(0, n_ests)),
                     meq = 1)
}
