combine_estimators <- function(ests, name_0 = NULL, boot_ests = NULL, cov = NULL, print = FALSE, exclude_t0 = FALSE, bias_type = 'raw_diff', ate_0 = NULL, ...) {
  # browser()
  if (is.null(boot_ests) & is.null(cov)) {
    stop("Must enter either resampled estimates or covariance estimate.")
  }
  if (!is.null(boot_ests)) {
    boot_ests[boot_ests == Inf] <- NA
    boot_ests[boot_ests == -Inf] <- NA
    use_i <- as.vector(!is.na(ests) & colMeans(is.na(boot_ests)) < 0.5)
    rm_b <- boot_ests[,use_i]

    C <- cov(rm_b, use = 'pairwise')
    mean_boot_ests <- apply(rm_b, 2, function(x) mean(x, na.rm = TRUE))
  }
  if (!is.null(cov)) {
    use_i <- as.vector(!is.na(ests) & colMeans(is.na(cov)) < 1)
    C <- cov[use_i,use_i]
  }

  C <- round(C, 10)
  if (!matrixcalc::is.positive.definite(round(C, 10))) {
    C <- Matrix::nearPD(C)$mat %>% as.matrix
  }
  

  rm_ests <- ests[,use_i]
  n_ests <- length(rm_ests)

  est_names <- colnames(rm_ests)
  if (is.null(name_0)) {
    synths <- lapply(est_names, function(n) {
      do_combination(ests = rm_ests, name_0 = n, C = C, 
                     exclude_t0 = exclude_t0, bias_type = bias_type,
                     boot_mean = mean_boot_ests,
                     ate_0 = ate_0)
      })
    synth_ates <- map(synths, 'synthetic_ate') %>% unlist
    shrunk_ates <- map(synths, 'shrunk_ate') %>% unlist
    synth_mse <- map(synths, 'naive_mse') %>% unlist
    shrunk_mse <- map(synths, 'shrunk_mse') %>% unlist
    synth_mse2 <- map(synths, 'naive_mse2') %>% unlist
    shrunk_mse2 <- map(synths, 'shrunk_mse2') %>% unlist
    # thr_var <- map(synths, 'th_var') %>% unlist

    synth_b <- map(synths, 'b') %>% simplify
    shrunk_b <- map(synths, 'b_shrink') %>% simplify
    shrinkage_facs <- map(synths, 'shrinkage_factor') %>% simplify

    ate_res <- bind_rows(
      list(
        data.frame(
          ate = synth_ates,
          theta_0 = est_names,
          synthetic = TRUE,
          var = synth_mse,
          # thr_var = thr_var,
          # var2 = synth_mse2,
          shrunk = FALSE
        ),
        data.frame(
          ate = shrunk_ates,
          theta_0 = est_names,
          synthetic = TRUE,
          var = shrunk_mse,
          # var2 = shrunk_mse2,
          shrunk = TRUE
        )
      ))
    # browser()
    en_mat <- matrix(est_names, length(est_names), length(est_names))
    if (exclude_t0) {
      use_names <- en_mat[row(en_mat) != col(en_mat)]
      ln <- length(est_names) - 1
    } else {
      use_names <- est_names
      ln <- length(est_names)
    }
    b_res <- bind_rows(
      list(
        data.frame(
          b = synth_b,
          est = use_names,
          theta_0 = rep(est_names, each = ln),
          shrunk = FALSE,
          shrinkage_factor = 0
        ),
        data.frame(
          b = shrunk_b,
          est = use_names,
          theta_0 = rep(est_names, each = ln),
          shrunk = TRUE,
          shrinkage_factor = shrinkage_facs
        )
      )
    )
  } else {
    comb <- do_combination(ests = rm_ests, name_0 = name_0, C = C, 
                           exclude_t0 = exclude_t0, bias_type = bias_type,
                           boot_mean = mean_boot_ests,
                           ate_0 = ate_0)
    synth_ates <- comb$synthetic_ate
    shrunk_ates <- comb$shrunk_ate
    synth_mse <- comb$naive_mse
    shrunk_mse <- comb$shrunk_mse
    synth_mse2 <- comb$naive_mse2
    shrunk_mse2 <- comb$shrunk_mse2

    synth_b <- comb$b
    shrunk_b <- comb$b_shrink
    shrinkage_facs <- comb$shrinkage_factor

    ate_res <- bind_rows(
      list(
        data.frame(
          ate = synth_ates,
          theta_0 = name_0,
          synthetic = TRUE,
          var = synth_mse,
          # var2 = synth_mse2,
          shrunk = FALSE
        ),
        data.frame(
          ate = shrunk_ates,
          theta_0 = name_0,
          synthetic = TRUE,
          var = shrunk_mse,
          # thr_var = thr_var,
          # var2 = shrunk_mse2,
          shrunk = TRUE
        )
      ))
    b_res <- bind_rows(
      list(
        data.frame(
          b = synth_b,
          est = est_names,
          theta_0 = name_0,
          shrunk = FALSE,
          shrinkage_factor = 0
        ),
        data.frame(
          b = shrunk_b,
          est = est_names,
          theta_0 = name_0,
          shrunk = TRUE,
          shrinkage_factor = shrinkage_facs
        )
      )
    )
  }




  list(ate_res = ate_res,
       b_res = b_res,
       C = C
       )
}


combine_estimators_linear <- function(data, fn_0, other_fns, print = FALSE, ...) {
  # browser()
  omega <- lapply(other_fns, function(f) {
    f(data, ...)
  }) %>% bind_cols
  n_ests <- length(other_fns) + 1
  omega_0 <- fn_0(data, ...) %>% unlist

  all_omegas <- cbind(omega_0, omega)


  # C_omega <- cov(all_omegas*data$y)
  C <- t(as.matrix(all_omegas)) %*% as.matrix(all_omegas) * var(data$y)
  B <- colSums((all_omegas - omega_0)*data$y)
  n <- nrow(data)

  qq <- C + B %*% t(B)

  Amat <- cbind(rep(1, n_ests), diag(n_ests))
  convex_soln <- quadprog::solve.QP(qq, dvec = rep(0, n_ests),
                                    Amat = Amat, bvec = c(1, rep(0, n_ests)),
                                    meq = 1)
  if (print) print(convex_soln)
  b_convex <- convex_soln$solution

  convex_omega <- c( as.matrix(all_omegas) %*% b_convex)
  dumb_omega <- all_omegas %>% rowMeans

  convex_ate <- sum(data$y*convex_omega)
  dumb_ate <-  sum(data$y*dumb_omega)
  list(convex_ate = convex_ate,
       dumb_ate = dumb_ate,
       C = C,
       # C_omega = C_omega,
       B = B,
       b = b_convex)
}

do_combination <- function(ests, name_0, C, print = FALSE, exclude_t0 = FALSE, is_cv = FALSE, bias_type = 'raw_diff', boot_mean = NULL, ate_0 = NULL) {
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
      B <- boot_mean - ate_0
    } else B <- boot_mean - ate_0
  } else if (bias_type == 'bootstrap_all') {
    if (is.null(boot_mean)) stop('Need bootstrap samples to compute bootstrap bias type.')
    B <- boot_mean - boot_mean[i_0]
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




  if (print) print(convex_soln)
  b_convex <- convex_soln$solution
  convex_ate <- unlist(ests) %*% b_convex
  b_shrink <- shrinkage_soln$solution
  shrunk_ate <- unlist(ests) %*% b_shrink

  r_mat <- matrix(r, n_ests-1, n_ests-1)
  pn <- v_0 - r
  tn <- v_0 + C[-i_0,-i_0] - r_mat - t(r_mat)
  w_aff <- solve(tn + B[-i_0] %*% t(B[-i_0])) %*% pn

  v_0 - t(pn) %*% solve(tn) %*% pn

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
  k <- n_ests-1
  omega <- rchisq(10000, df = k)
  if (is_cv) {
    qty <- mean(omega^2/k/(1+omega)^2)
  } else qty <- mean(omega^3/k/(1+omega)^2)
  var_est <- v_0 - t(pn) %*% solve(tn) %*% pn *(1 - qty)

  # browser()
  naive_var <- t(b_convex) %*% C %*% b_convex
  naive_mse <- naive_var + (t(b_convex) %*% B)^2
  naive_mse2 <- (sqrt(naive_var) + 1/2*abs(t(b_convex) %*% B))^2

  shrunk_var <- t(b_shrink) %*% C %*% b_shrink
  shrunk_mse <- shrunk_var + (t(b_shrink) %*% B)^2
  shrunk_mse2 <- (sqrt(shrunk_var) + 1/2*abs(t(b_shrink) %*% B))^2
  list(b = b_convex, synthetic_ate = convex_ate, b_shrink = b_shrink,
       shrunk_ate = shrunk_ate, shrinkage_factor = w,
       naive_var = naive_var, naive_mse = naive_mse,
       naive_mse2 = naive_mse2,
       shrunk_var = shrunk_var, shrunk_mse = shrunk_mse,
       shrunk_mse2 = shrunk_mse2,
       th_var = var_est)
}

qp <- function(qq, n_ests) {
  Amat <- cbind(rep(1, n_ests), diag(n_ests))
  quadprog::solve.QP(qq, dvec = rep(0, n_ests),
                                    Amat = Amat, bvec = c(1, rep(0, n_ests)),
                                    meq = 1)
}

get_deriv <- function(all_theta, v_0, c, r, B, b, i_0) {
  # browser()
  theta <- all_theta[-i_0]
  theta_0 <- all_theta[i_0]
  p <- length(r)
  bone <- rep(1, p)
  sigmatilde <- matrix(v_0, p, p) - matrix(r, p, p) - matrix(r, p, p, byrow = TRUE) - c
  num <- v_0 * bone - r
  bb <- B[B != 0] %*% t(B[B != 0])

  dB <- matrix(0, p+1, p+1)

  D_0 <- matrix(0, p, p)
  for (i in 1:p) {
    for (k in 1:p) {
      D_0[i,k] <- unlist(2*theta_0 - theta[i] - theta[k])
    }
  }
# browser()
  upper_right <- solve(sigmatilde + bb, D_0) %*% solve(sigmatilde + bb, num)
  upper_left <- -sum(upper_right)
  dB[i_0,1] <- upper_left
  dB[i_0,-1] <- upper_right

  for (j in 1:p) {
    D <- matrix(0, p, p)
    for (i in 1:p) {
      for (k in 1:p) {
        if (i == j & k == j) {
          D[i,k] <- unlist(2*(theta[j] - theta_0))
        } else if (i == j) {
          D[i,k] <- unlist(theta[k] - theta_0)
        } else if (k == j) {
          D[i,k] <- unlist(theta[i] - theta_0)
        }
      }
    }
    if (j < i_0) {
      dB[j,-1] <- solve(sigmatilde + bb, D) %*% solve(sigmatilde + bb, num)
      dB[j, 1] <- -sum(dB[j + 1,-1])
    } else {
      dB[j + 1,-1] <- solve(sigmatilde + bb, D) %*% solve(sigmatilde + bb, num)
      dB[j+1, 1] <- -sum(dB[j + 1,-1])
    }

  }
# browser()

  dB %*% unlist(all_theta) + b
}
