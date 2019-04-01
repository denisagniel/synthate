fit_fn <- function(ds, gen_mod, ate_list, B, C = NULL) {
  outcome_fm <- gen_mod$outcome_fm
  outcome_fam <- gen_mod$outcome_fam
  ps_fm <- gen_mod$ps_fm
  ps_fam <- gen_mod$ps_fam
  cov_ids <- gen_mod$cov_ids
  
  ds <- estimate_scores(ds, outcome_fm = outcome_fm,
                        ps_fm = ps_fm,
                        ps_fam = ps_fam,
                        outcome_fam = outcome_fam)
  
  print('scores estimated...')
  thetahat <- estimate_ates(ds,
                            ate_list,
                            cov_ids = cov_ids,
                            outcome_fm = stringr::str_c('d + ', outcome_fm),
                            outcome_fam = outcome_fam)
  print('initial thetas estimated...')
  X <- ds %>% select(one_of(cov_ids))
  W <- ds %>% pull(d)
  Y <- ds$y
  
  rf_fit <- randomForest::randomForest(as.formula(stringr::str_c('y ~ d + ', outcome_fm)), data = ds)
  grf_fit <- grf::causal_forest(X = X, Y = Y, W = W)
  # browser()
  print('prediction model fit...')
  # browser()
  predict_delta <- function(d) {
    as.vector(predict(grf_fit, newdata = d)$predictions)
  }
  predict_y <- function(d) {
    unlist(predict(rf_fit, newdata = d))
  }
  # browser()
  resample_thetas <- 
    resample_fn(dat = ds,
                ypredfn = predict_y,
                dpredfn = predict_delta,
                B = B,
                ate_list = ate_list,
                outcome_fm = outcome_fm,
                ps_fm = ps_fm,
                ps_fam = ps_fam,
                outcome_fam = outcome_fam)
  print('resampling done...')
  boot_theta <- resample_thetas[[1]] %>% as.matrix
  null_theta <- resample_thetas[[2]] %>% as.matrix
  
  
  this_cov <- cov(boot_theta)
  #------------------------
  ## this way uses regular bootstrap and raw differences
  #--------------------------
  old_way <- combine_estimators(thetahat,
                                boot_ests = boot_theta,
                                name_0 = 'ate_dr',
                                bias_type = 'raw_diff',
                                cov = C)
  
  #---------------------------
  ## de-meaning method
  null_way <- combine_estimators(thetahat,
                                 boot_ests = null_theta,
                                 name_0 = 'ate_dr',
                                 bias_type = 'bootstrap',
                                 ate_0 = 0,
                                 cov = C)
  
  #--------------------------
  ## just using the prediction model for the ate, but not for the bootstrapping
  gn_ate_0 <- mean(
    predict_delta(ds)
  )
  hybrid_boot_gn_way <- combine_estimators(thetahat,
                                           boot_ests = boot_theta,
                                           name_0 = 'ate_dr',
                                           bias_type = 'bootstrap',
                                           ate_0 = gn_ate_0,
                                           cov = C)
  #----------------------------
  ## shrunk bias
  #---------------------------
  shrunk_bias_way <- combine_estimators(thetahat,
                                        boot_ests = boot_theta,
                                        name_0 = 'ate_dr',
                                        bias_type = 'shrunk',
                                        n = nrow(ds), cov = C)
  
  alvvays <- list(
    old_way,
    null_way,
    hybrid_boot_gn_way,
    shrunk_bias_way
  )
  
  all_synthetic_thetas <- map(alvvays, 'ate_res') %>%
    bind_rows %>%
    filter(!shrunk) %>%
    mutate(type = c('old', 'null', 'hybrid_gn', 'shrunk'),
           true_ate = gen_mod$true_ate)
  
  all_regular_thetas <- gather(thetahat, theta_0, ate) %>%
    mutate(synthetic = FALSE,
           var = NA,
           shrunk = FALSE,
           type = 'raw',
           true_ate = gen_mod$true_ate)
  
  all_thetas <- all_regular_thetas %>%
    full_join(all_synthetic_thetas)
  all_bs <- map(alvvays, 'b_res') %>%
    bind_rows %>%
    filter(!shrunk) %>%
    mutate(type = rep(c('old', 'null', 'hybrid_gn', 'shrunk'), each = 4))
  
  return(list(thetas = all_thetas,
              bs = all_bs,
              covhat = this_cov))
}

split_sample_fn <- function(ds, gen_mod, ate_list, B, C = NULL) {
  a_data <- ds %>%
    sample_frac(0.5)
  b_data <- ds %>%
    anti_join(a_data)
  c(a_thetahat, a_bhat, a_Chat) %<-% fit_fn(a_data, 
                                            gen_mod = gen_mod, 
                                            ate_list = ate_list,
                                            C = C,
                                            B = B)
  c(b_thetahat, b_bhat, b_Chat) %<-% fit_fn(b_data, 
                                            gen_mod = gen_mod, 
                                            ate_list = ate_list,
                                            C = C,
                                            B = B)
  ab_ds <- a_thetahat %>%
    filter(type == 'raw') %>%
    select(theta_0, ate) %>%
    merge(b_bhat %>% select(-theta_0), by.x = 'theta_0', by.y = 'est') %>%
    group_by(type) %>%
    summarise(theta_s = sum(ate*b))
  
  ba_ds <- b_thetahat %>%
    filter(type == 'raw') %>%
    select(theta_0, ate) %>%
    merge(a_bhat %>% select(-theta_0), by.x = 'theta_0', by.y = 'est') %>%
    group_by(type) %>%
    summarise(theta_s = sum(ate*b))
  
  split_sample_ests <- inner_join(ab_ds, ba_ds, by = 'type') %>%
    transmute(type, theta_s.x, theta_s.y, theta_s = theta_s.x/2 + theta_s.y/2)
  
  split_sample_bs <- b_bhat %>%
    inner_join(a_bhat, by = c('est', 'theta_0', 'type')) %>%
    transmute(est, theta_0, type, b.x, b.y, b = b.x/2 + b.y/2)
  
  list(ss_thetas = split_sample_ests,
       ss_bs = split_sample_bs)
}
