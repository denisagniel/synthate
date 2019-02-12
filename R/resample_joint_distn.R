resample_joint_distn <- function(dat, ypredfn, dpredfn) {
  mod_dat <- dat %>%
    mutate(yhat = ypredfn(dat),
           e = y - yhat,
           ate_0 = dpredfn(dat))
  e <- mod_dat$e
  mod_dat <- mod_dat %>% select(-e)
  new_dat <- mod_dat %>%
    sample_frac(replace = TRUE) %>%
    mutate(e = sample(e, replace = TRUE),
           y_e = yhat + e)
  new_dat
}

resample_fn <- function(dat, ypredfn, dpredfn, B,
                        ate_list,
                        outcome_fm,
                        ps_fm,
                        ps_fam,
                        outcome_fam) {
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
        d == 1 ~ y - ate_0
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
  return(map(list(theta_boot = theta_boot, 
                  theta_adj = theta_adj, 
                  theta_bing = theta_bing), bind_rows))
}

