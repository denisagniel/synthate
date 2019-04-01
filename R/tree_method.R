tree_method_bhat <- function(dat, theta_0, B = 200, ate_list, cov_ids, outcome_fm, outcome_fam) {
  mod_data <- remove_effect_from_data(dat, theta_0)
  n <- nrow(mod_data)
  #----------------------------------
  ## bootstrap mod data
  #----------------------------------
  boot_l <- list()
  for (b in 1:B) {
    boot_data <- mod_data %>%
      sample_n(n, replace = TRUE)
    boot_theta <- estimate_ates(boot_data,
                                ate_list,
                                cov_ids = cov_ids,
                                outcome_fm = stringr::str_c('d + ', outcome_fm),
                                outcome_fam = outcome_fam)
    boot_l[[b]] <- boot_theta
  }
  boot_ests <- bind_rows(boot_l) %>% as.matrix
  mse_b <- t(boot_ests) %*% boot_ests
  convex_b <- qp(mse_b, nrow(mse_b))$solution
  convex_b
}
