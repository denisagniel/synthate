resample_joint_distn <- function(dat, predfn) {
  yhat <- predfn(dat)
  mod_dat <- dat %>%
    mutate(yhat = predfn(dat),
           e = y - yhat,
           y_0 = predfn(dat %>% mutate(d = 0)),
           y_1 = predfn(dat %>% mutate(d = 1)),
           ate_0 = mean(y_1 - y_0))
  e <- mod_dat$e
  mod_dat <- mod_dat %>% select(-e)
  new_dat <- mod_dat %>%
    sample_frac(replace = TRUE) %>%
    mutate(e = sample(e, replace = TRUE),
           y_e = yhat + e)
  new_dat
}

resample_fn <- function(dat, predfn, B) {
  theta_bing <- theta_adj <- theta_boot <- list()
  for (b in 1:B) {
    resample_data <- resample_joint_distn(dat, predfn)
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
    
    resample_adj <- resample_data %>%
      mutate(y = case_when(
        d == 0 ~ y,
        d == 1 ~ y - y_1 + y_0
      ))
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
  return(list(theta_boot, theta_adj, theta_bing))
}

