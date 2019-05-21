#' Resample ATEs.
#'
#' @param dat data frame for analysis.
#' @param ypredfn function to predict y, mostly deprecated, default is NULL.
#' @param dpredfn function to predict the ATE used to demean data.
#' @param B number of bootstrap replications.
#' @param ate_list list of ATE functions.
#' @param outcome_fm outcome formula.
#' @param ps_fm propensity score formula.
#' @param ps_fam propensity score family.
#' @param outcome_fam outcome family.
#' @param cov_ids names of covariates for matching/balancing.
#'
#' @return list of tibbles of resampled data.
#' @export
#'
#' @examples
#' @importFrom stringr str_c
resample_fn <- function(dat, ypredfn = NULL, dpredfn, B,
                        ate_list,
                        outcome_fm,
                        ps_fm,
                        ps_fam,
                        outcome_fam,
                        cov_ids) {
  theta_bing <- theta_adj <- theta_boot <- list()
  for (b in 1:B) {
    resample_data <- resample_joint_distn(dat, ypredfn, dpredfn)
    resample_data <- estimate_scores(resample_data, 
                                     outcome_fm = outcome_fm,
                                     ps_fm = ps_fm,
                                     ps_fam = ps_fam,
                                     outcome_fam = outcome_fam)
    theta_boot[[b]] <- estimate_ates(resample_data,
                                     ate_list,
                                     cov_ids = cov_ids,
                                     outcome_fm = stringr::str_c('d + ', 
                                                                 outcome_fm),
                                     outcome_fam = outcome_fam)
    
    if (!is.null(ypredfn)) {
      resample_bing <- resample_data %>%
        mutate(y = y_e)
      resample_bing <- estimate_scores(resample_bing, 
                                       outcome_fm = outcome_fm,
                                       ps_fm = ps_fm,
                                       ps_fam = ps_fam,
                                       outcome_fam = outcome_fam)
      theta_bing[[b]] <- estimate_ates(resample_bing,
                                       ate_list,
                                       cov_ids = cov_ids,
                                       outcome_fm = stringr::str_c('d + ', 
                                                                   outcome_fm),
                                       outcome_fam = outcome_fam)
      
    }
    # browser()
    #---------------------
    ## horrible kluge to account for waernbaum
    #---------------------
    if (is.integer(resample_data$y)) {
      resample_adj <- resample_data %>%
        mutate(y = case_when(
          d == 0 ~ as.double(y),
          d == 1 ~ y - ate_0
        ),
        y = if_else(y < 0, 0, y))
    } else {
      resample_adj <- resample_data %>%
        mutate(y = case_when(
          d == 0 ~ y,
          d == 1 ~ y - ate_0
        ))
    }
    
    
    resample_adj <- estimate_scores(resample_adj, 
                                    outcome_fm = outcome_fm,
                                    ps_fm = ps_fm,
                                    ps_fam = ps_fam,
                                    outcome_fam = outcome_fam)
    theta_adj[[b]] <- estimate_ates(resample_adj,
                                    ate_list,
                                    cov_ids = cov_ids,
                                    outcome_fm = stringr::str_c('d + ', 
                                                                outcome_fm),
                                    outcome_fam = outcome_fam)
  }
  if (!is.null(ypredfn)) {
    return(map(list(theta_boot = theta_boot, 
                    theta_adj = theta_adj, 
                    theta_bing = theta_bing), bind_rows))
  } else {
    return(map(list(theta_boot = theta_boot, 
                    theta_adj = theta_adj), bind_rows))
  }
}

resample_joint_distn <- function(dat, ypredfn = NULL, dpredfn) {
  mod_dat <- dat %>%
    mutate(ate_0 = dpredfn(dat))
  if (!is.null(ypredfn)) {
    mod_dat %>%
      mutate(yhat = ypredfn(dat),
             e = y - yhat)
    e <- mod_dat$e
    mod_dat <- mod_dat %>% select(-e)
    new_dat <- mod_dat %>%
      sample_frac(replace = TRUE) %>%
      mutate(e = sample(e, replace = TRUE),
             y_e = yhat + e)
  } else{
    new_dat <- mod_dat %>%
      sample_frac(replace = TRUE)
  }
  new_dat
}



