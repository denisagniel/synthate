#' Workhorse function for creating a synthetic estimator by combining multiple candidate estimators.
#'
#' Creates a synthetic estimator by minimizing the (estimated) mean squared error of a linear combination of multiple candidate estimators.
#'
#' @param ests one-row, p-column data frame of estimators.
#' @param name_0 character value of the name of the presumed unbiased estimator, \eqn{\theta_0}. Default is NULL, which returns results for using each candidate estimator as \eqn{\theta_0}, one synthetic estimator for each. If \code{ate_0} is given, then it is used as \eqn{\theta_0} in place of this.
#' @param C p x p covariance matrix of \code{ests}.
#' @param print logical indicating whether details should be printed. Default is FALSE.
#' @param exclude_t0 logical indicating whether \eqn{\theta_0} should be considered an external estimator (not a candidate for combining with others). Default is FALSE.
#' @param bias_type method to compute the bias in the mean squard error. Default is \code{raw_diff}, which computes the bias as the raw difference between each of the candidate estimator and \eqn{\theta_0}. Other options to compute the bias include: \code{bootstrap} which computes the bias as the difference between the mean of the bootstrap samples and the observed value of \eqn{\theta_0}; \code{bootstrap_all} which computes the bias as the mean of the difference between the bootstrapped version of the candidate estimator and the bootstrapped version of \eqn{\theta_0}; \code{none} which assumes no bias; \code{shrunk} which computes the bias as the raw difference divided by \code{n}. 
#' @param boot_mean mean of bootstrap samples for bootstrap-based bias estimation.
#' @param ate_0 external value of \eqn{\theta_0}. Default is NULL, in which case \eqn{\theta_0} is taken to be \code{name_0}. 
#' @param n sample size. Default is NULL. Needed only if \code{bias_type} is \code{shrunk}.
#'
#' @return list of three objects, including \code{ate_res} which gives results for the synthetic estimator, \code{b_res} which gives results for how the estimators were combined, and \code{C} which gives the covariance matrix of the estimators. 
#' 
#' @import dplyr

do_combination <- function(ests, name_0, C, print = FALSE, exclude_t0 = FALSE, bias_type = 'raw_diff', boot_mean = NULL, ate_0 = NULL, n = NULL) {
  # browser()
  est_0 <- ests %>% select_(name_0) %>% unlist
  i_0 <- which(colnames(ests) == name_0)
  v_0 <- C[i_0, i_0]
  v <- diag(C)[-i_0]
  r <- C[i_0,][-i_0]
  
  if (exclude_t0) {
    est_names <- colnames(ests)
    est_names <- est_names[est_names != name_0]
    ests <- ests %>% select_(.dots = est_names)
    C <- C[-i_0, -i_0]
  }
  n_ests <- length(ests)
  
  if (bias_type == 'raw_diff') {
    #-------------------------
    # raw differences
    if (is.null(ate_0)) {
      B <- unlist(ests) - unlist(est_0)
    } else B <- unlist(ests) - ate_0
    
  } else if (bias_type == 'bootstrap') {
    if (is.null(boot_mean)) stop('Need bootstrap samples to compute bootstrap bias type.')
    
    if (is.null(ate_0)) {
      B <- boot_mean - unlist(est_0)
    } else B <- boot_mean - ate_0
  } else if (bias_type == 'bootstrap_all') {
    if (is.null(boot_mean)) stop('Need bootstrap samples to compute bootstrap bias type.')
    B <- boot_mean - boot_mean[i_0]
  } else if (bias_type == 'none') {
    B <- rep(0, n_ests)
  } else if (bias_type == 'shrunk') {
    if (is.null(n)) stop('n must be provided for shrunk bias.')
    if (is.null(ate_0)) {
      B <- (unlist(ests) - unlist(est_0))/n
    } else B <- (unlist(ests) - ate_0)/n
  }
  
  qq <- C + B %*% t(B)
  qq_adj <- qq/norm(qq, '2')
  convex_soln <- qp(qq_adj, n_ests)
  # cs_adj <- qp(qq_adj, n_ests)
  
  #-------------------------
  # adjusted differences
  w <- B
  w[-i_0] <- B[-i_0]^2/(v_0 + v - 2*r + B[-i_0]^2)
  Btilde <- B*w
  qqtilde <- C + Btilde %*% t(Btilde)
  qqt_adj <- qqtilde/norm(qqtilde, '2')
  shrinkage_soln <- qp(qqt_adj, n_ests)
  
  msehat <- asymptotic_mse(V_0 = v_0,
                           V = C[-i_0,-i_0],
                           C = r,
                           deltahat = B[-i_0])
  shrunk_msehat <- asymptotic_mse(V_0 = v_0,
  V = C[-i_0,-i_0],
  C = r,
  deltahat = Btilde[-i_0])
  
  
  if (print) print(convex_soln)
  b_convex <- convex_soln$solution
  convex_ate <- unlist(ests) %*% b_convex
  b_shrink <- shrinkage_soln$solution
  shrunk_ate <- unlist(ests) %*% b_shrink
  
  # r_mat <- matrix(r, n_ests-1, n_ests-1)
  # pn <- v_0 - r
  # tn <- v_0 + C[-i_0,-i_0] - r_mat - t(r_mat)
  # w_aff <- solve(tn + B[-i_0] %*% t(B[-i_0])) %*% pn
  
  # v_0 - t(pn) %*% solve(tn) %*% pn
  
  # gamma <- mvnfast::rmvn(1000, mu = rep(0, n_ests-1), sigma = tn)
  # g <- diag(gamma %*% solve(tn) %*% t(gamma))
  # qty <- rep(0, 1000)
  # for (i in 1:1000) {
  #   gamma_i <- gamma[i,]
  #   g_i <- t(gamma_i) %*% solve(tn) %*% gamma_i
  #   qty[i] <- g_i^2/(1+g_i)^2 * t(pn) %*% solve(tn) %*% gamma_i %*% t(gamma_i) %*% solve(tn) %*% pn
  # }
  # browser()
  # var_est <- v_0 - t(pn) %*% solve(tn) %*% pn + mean(qty)
  # k <- n_ests-1
  # omega <- rchisq(10000, df = k)
  # if (is_cv) {
  #   qty <- mean(omega^2/k/(1+omega)^2)
  # } else qty <- mean(omega^3/k/(1+omega)^2)
  # var_est <- v_0 - t(pn) %*% solve(tn) %*% pn *(1 - qty)
  
  # browser()
  # naive_var <- t(b_convex) %*% C %*% b_convex
  # naive_mse <- naive_var + (t(b_convex) %*% B)^2
  # naive_mse2 <- (sqrt(naive_var) + 1/2*abs(t(b_convex) %*% B))^2
  # 
  # shrunk_var <- t(b_shrink) %*% C %*% b_shrink
  # shrunk_mse <- shrunk_var + (t(b_shrink) %*% B)^2
  # shrunk_mse2 <- (sqrt(shrunk_var) + 1/2*abs(t(b_shrink) %*% B))^2
  list(b = b_convex, 
       synthetic_ate = convex_ate, 
       b_shrink = b_shrink,
       shrunk_ate = shrunk_ate,
       # shrinkage_factor = w,
       # naive_var = naive_var, 
       # naive_mse = naive_mse,
       # naive_mse2 = naive_mse2,
       asymp_mse = msehat,
       shrunk_mse = shrunk_msehat
       # shrunk_var = shrunk_var, 
       # shrunk_mse = shrunk_mse,
       # shrunk_mse2 = shrunk_mse2,
       # th_var = var_est
       )
}
