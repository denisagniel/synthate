#' Create a synthetic estimator by combining multiple candidate estimators.
#'
#' Creates a synthetic estimator by minimizing the (estimated) mean squared error of a linear combination of multiple candidate estimators.
#'
#' @param ests one-row, p-column data frame of estimators.
#' @param name_0 character value of the name of the presumed unbiased estimator, \eqn{\theta_0}. Default is NULL, which returns results for using each candidate estimator as \eqn{\theta_0}, one synthetic estimator for each. If \code{ate_0} is given, then it is used as \eqn{\theta_0} in place of this.
#' @param boot_ests p-column matrix of bootstrap estimators corresponding to the estimators in \code{ests}. If NULL, then \code{cov} must be supplied.
#' @param cov p x p covariance matrix of \code{ests}. Default is NULL, but if supplied \code{boot_ests} are not needed.
#' @param print logical indicating whether details should be printed. Default is FALSE.
#' @param exclude_t0 logical indicating whether \eqn{\theta_0} should be considered an external estimator (not a candidate for combining with others). Default is FALSE.
#' @param bias_type method to compute the bias in the mean squard error. Default is \code{raw_diff}, which computes the bias as the raw difference between each of the candidate estimator and \eqn{\theta_0}. Other options to compute the bias include: \code{bootstrap} which computes the bias as the difference between the mean of the bootstrap samples and the observed value of \eqn{\theta_0}; \code{bootstrap_all} which computes the bias as the mean of the difference between the bootstrapped version of the candidate estimator and the bootstrapped version of \eqn{\theta_0}; \code{none} which assumes no bias; \code{shrunk} which computes the bias as the raw difference divided by \code{n}. 
#' @param ate_0 external value of \eqn{\theta_0}. Default is NULL, in which case \eqn{\theta_0} is taken to be \code{name_0}. 
#' @param n sample size. Default is NULL. Needed only if \code{bias_type} is \code{shrunk}.
#'
#' @return list of three objects, including \code{ate_res} which gives results for the synthetic estimator, \code{b_res} which gives results for how the estimators were combined, and \code{C} which gives the covariance matrix of the estimators. 
#'
#' @examples
#'
#' gen_mod <- generate_data(n = 100, 
#'                          dgp = 'ks', 
#'                          correct_outcome = FALSE,
#'                          correct_ps = TRUE)
#' this_data <- gen_mod$data
#' this_data <- estimate_scores(this_data, outcome_fm = outcome_fm,
#'                              ps_fm = ps_fm,
#'                              ps_fam = ps_fam,
#'                              outcome_fam = outcome_fam)
#' thetahat <- with(gen_mod,
#'                  estimate_ates(this_data,
#'                           ate_list,
#'                           cov_ids = cov_ids,
#'                           outcome_fm = stringr::str_c('d + ', outcome_fm),
#'                           outcome_fam = outcome_fam))
#'                           predict_delta <- function(d) {
#'                           gen_mod$true_ate
#'                           }
#'  resample_thetas <- with(gen_mod, resample_fn(dat = this_data,
#'                                 dpredfn = predict_delta,
#'                                 B = B,
#'                                 ate_list = ate_list,
#'                                 outcome_fm = outcome_fm,
#'                                 ps_fm = ps_fm,
#'                                 ps_fam = ps_fam,
#'                                 outcome_fam = outcome_fam))
#'  boot_theta <- resample_thetas[[1]] %>% as.matrix
#'  
#'  synthetic_estimator <- combine_estimators(thetahat,
#'  boot_ests = boot_theta,
#'  name_0 = 'ate_dr',
#'  bias_type = 'raw_diff')
#'  synthetic_estimator$ate_res
#' @export 

combine_estimators <- function(ests, 
                               name_0 = NULL, 
                               boot_ests = NULL, 
                               cov = NULL, 
                               print = FALSE, 
                               exclude_t0 = FALSE, 
                               bias_type = 'raw_diff', 
                               ate_0 = NULL, 
                               n = NULL,...) {
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
                     ate_0 = ate_0,
                     n = n)
      })
    synth_ates <- map(synths, 'synthetic_ate') %>% unlist
    shrunk_ates <- map(synths, 'shrunk_ate') %>% unlist
    synth_mse <- map(synths, 'naive_mse') %>% unlist
    shrunk_mse <- map(synths, 'shrunk_mse') %>% unlist
    synth_mse2 <- map(synths, 'naive_mse2') %>% unlist
    asymp_mse <- map(synths, 'asymp_mse') %>% unlist
    # shrunk_mse2 <- map(synths, 'shrunk_mse2') %>% unlist
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
          var = asymp_mse,
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
                           ate_0 = ate_0, n = n)
    synth_ates <- comb$synthetic_ate
    shrunk_ates <- comb$shrunk_ate
    synth_mse <- comb$naive_mse
    shrunk_mse <- comb$shrunk_mse
    synth_mse2 <- comb$naive_mse2
    shrunk_mse2 <- comb$shrunk_mse2
    asymp_mse <- comb$asymp_mse

    synth_b <- comb$b
    shrunk_b <- comb$b_shrink
    shrinkage_facs <- comb$shrinkage_factor

    ate_res <- bind_rows(
      list(
        data.frame(
          ate = synth_ates,
          theta_0 = name_0,
          synthetic = TRUE,
          var = asymp_mse,
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

