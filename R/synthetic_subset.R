#' Compute synthetic estimator on a subset of candidate estimators.
#' 
#' This function allows you to compute a set of candidate estimators and bootstrap estimates once and then compute on it many times. It's designed to work well with \code{map2}.
#'
#' @param theta one-row data frame of candidate estimators.
#' @param boot data frame of bootstrapped versions of candidate estimators.
#' @param estimators vector of estimator names to be included in synthetic estimator; if NULL (the default) all estimators are used.
#' @param ... additional arguments to be passed to \code{combine_estimators}.
#'
#' @return vector of synthetic estimators
#' @export
#'
#' @import dplyr
#' @importFrom purrr simplify
#'
#' @examples
#' 
#' ate_list <- list(
#' ipw1 = ipw1_ate,
#' ipw2 = ipw2_ate,
#' ipw3 = ipw3_ate,
#' strat = strat_ate,
#' strat_regr = strat_regr_ate,
#' match_ps = match_ps_ate,
#' match_prog = match_prog_ate,
#' match_both = match_both_ate,
#' cal_ps = caliper_ps_ate,
#' cal_both = caliper_both_ate,
#' regr = regr_ate,
#' dr = dr_ate,
#' bal = bal_ate,
#' hdbal = highdim_bal_ate
#' )
#' gen_mod <- generate_data(n = 100, 
#'                          dgp = 'ks', 
#'                          correct_outcome = FALSE,
#'                          correct_ps = TRUE)
#' this_data <- gen_mod$data
#' this_data <- with(gen_mod, estimate_scores(this_data, outcome_fm = outcome_fm,
#'                              ps_fm = ps_fm,
#'                              ps_fam = ps_fam,
#'                              outcome_fam = outcome_fam))
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
#'                                 B = 50,
#'                                 ate_list = ate_list,
#'                                 outcome_fm = outcome_fm,
#'                                 ps_fm = ps_fm,
#'                                 ps_fam = ps_fam,
#'                                 outcome_fam = outcome_fam,
#'                                 cov_ids = cov_ids))
#'  boot_theta <- resample_thetas[[1]]
#'  
#'  synthetic_subset(thetahat, boot_theta,
#'                   estimators = c('ate_regr', 'ate_ipw_2', 'ate_dr', 'ate_bal'),
#'                   name_0 = 'ate_dr')
#'  synthetic_subset(thetahat, boot_theta,
#'                   estimators = c('ate_regr', 'ate_ipw_2', 'ate_dr', 'match_ps'),
#'                   ate_0 = thetahat$ate_bal,
#'                   name_0 = 'ate_dr')
#'  
#'  out_df <- tibble::tibble(
#'  thetahat = list(thetahat),
#'  boot_theta = list(boot_theta))
#'  
#'  
#'  ## 
#'  
#'  dplyr::mutate(out_df, 
#'  theta_s = unlist(purrr::map2(
#'     thetahat, 
#'     boot_theta, 
#'     synthetic_subset, 
#'     estimators = c('ate_regr', 'ate_ipw_2', 'ate_dr', 'ate_bal'),
#'     name_0 = 'ate_dr')))
#'     
#'     
synthetic_subset <- function(theta, boot, estimators = NULL, ...) {
  if (!is.null(estimators)) {
    theta_use <- select(theta, one_of(estimators))
    boot_use <- select(boot_theta, one_of(estimators))
    boot_use <- as.matrix(boot_use)
  } else {
    theta_use <- theta
    boot_use <- as.matrix(boot)
  }
  # browser()
  comb <- combine_estimators(ests = theta_use,
                             boot_ests = boot_use,
                             ...)
  filter(comb$ate_res, !shrunk)
}
