psrg_sim_fn <- function(n,
                   p,
                   delta = 0.5,
                   compliance_p = 0.75,
                   fn_list, 
                   psfm,
                   run,
                   tmpdir) {
  # browser()
  library(purrr)
  library(tidyr)
  library(dplyr)
  library(clustermq)
  library(synthate)
  library(glue)
  
  bootfn <- function(x, i) {
    # browser()
    new_data <- x[i,]
    
    ps_model <- glm(psfm, family = binomial, data = new_data %>% filter(z == 1))
    new_data <- new_data %>%
      mutate(pr_score = predict(ps_model, newdata = new_data, type = 'response'),
             ps_grp = Hmisc::cut2(pr_score, g = 5),
             pi_c = mean(pr_score)) %>%
      group_by(ps_grp) %>%
      mutate(pi_cj = mean(pr_score)) %>%
      ungroup
    
    
    out <- estimate_ates(new_data,  fn_list)
    as.matrix(out)
  }
  
  set.seed(run)
  print(glue('Simulation {run} for n = {n}, p = {p}, compliance_p = {compliance_p}, delta = {delta}.'))
  sim_data <- 
    generate_data_psrg(n, p, delta, compliance_p)
  train_data <- sim_data %>% sample_frac(0.5)
  test_data <- sim_data %>% anti_join(train_data)
  
  #------------------------------
  ## first, all analysis on full data
  #-----------------------------
  ps_model <- glm(psfm, family = binomial, data = sim_data %>% filter(z == 1))
  
  sim_data <- sim_data %>%
    mutate(pr_score = predict(ps_model, newdata = sim_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  full_caces <- estimate_ates(sim_data, fn_list)
  b_theta <- boot::boot(sim_data, bootfn, R = 200)$t
  boot_theta_s <- combine_estimators(ests = full_caces, boot_ests = b_theta)
  full_synth_ests <- boot_theta_s$ate_res %>%
    group_by(theta_0, shrunk, synthetic) %>%
    transmute(estimate = ate,
              var = var,
              sample = 'full') %>%
    ungroup
  
  #----------------------------
  ## now use train data to fit model and combine on test data
  #---------------------------
  ps_model_1 <- glm(psfm, family = binomial, data = train_data %>% filter(z == 1))
  train_data <- train_data %>%
    mutate(pr_score = predict(ps_model_1, newdata = train_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  test_data <- test_data %>%
    mutate(pr_score = predict(ps_model_1, newdata = test_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  train_caces_1 <- estimate_ates(train_data, fn_list)
  test_caces_1 <- estimate_ates(test_data, fn_list)
  
  # b_theta_1 <- boot::boot(train_data, bootfn, R = 200)$t
  boot_theta_s_1 <- combine_estimators(ests = train_caces_1, boot_ests = b_theta)
  
  
  # bhat_1 <- boot_theta_s_1$b_res %>% filter(!shrunk) %>% select(b) %>% unlist
  # holdout_est_1 <- unlist(test_caces_1) %*% bhat_1
  
  test_caces_ds <- test_caces_1 %>%
    gather(est, cace)
  holdout_ests_1 <- test_caces_ds %>%
    inner_join(boot_theta_s_1$b_res) %>%
    group_by(theta_0, shrunk) %>%
    summarise(estimate = sum(cace*b))%>%
    mutate(sample = 'CV',
           synthetic = TRUE)
  
  #-------------------------------
  ## flip the script
  #-----------------------------
  ps_model_2 <- glm(psfm, family = binomial, data = test_data %>% filter(z == 1))
  train_data <- train_data %>%
    mutate(pr_score = predict(ps_model_2, newdata = train_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  test_data <- test_data %>%
    mutate(pr_score = predict(ps_model_2, newdata = test_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  train_caces_2 <- estimate_ates(train_data, fn_list)
  test_caces_2 <- estimate_ates(test_data, fn_list)
  
  # b_theta_2 <- boot::boot(test_data, bootfn, R = 200)$t
  boot_theta_s_2 <- combine_estimators(ests = test_caces_2, boot_ests = b_theta)
  
  # bhat_2 <- boot_theta_s_2$b_res %>% filter(!shrunk) %>% select(b) %>% unlist
  # holdout_est_2 <- unlist(train_caces_2) %*% bhat_2
  train_caces_ds <- train_caces_2 %>%
    gather(est, cace)
  holdout_ests_2 <- train_caces_ds %>%
    inner_join(boot_theta_s_2$b_res) %>%
    group_by(theta_0, shrunk) %>%
    summarise(estimate = sum(cace*b)) %>%
    mutate(sample = 'CV',
           synthetic = TRUE) %>%
    ungroup
  
  holdout_ests <- holdout_ests_1 %>%
    full_join(holdout_ests_2) %>%
    group_by(theta_0, shrunk, sample, synthetic) %>%
    summarise(estimate = mean(estimate)) %>%
    ungroup
  
  long_caces  <- full_caces %>%
    gather(theta_0, estimate) %>%
    mutate(var = diag(cov(b_theta)),
           shrunk = FALSE,
           sample = 'full',
           synthetic = FALSE) %>%
    full_join(full_synth_ests) %>%
    full_join(holdout_ests)
  
  out <- long_caces %>%
    mutate(n = n,
           p = p,
           compliance_p = compliance_p,
           delta = delta)
  saveRDS(out, 
          glue('{tmpdir}psrg-compliance-res-n{n}-p{p}-delta{delta}-compliance_p{compliance_p}-sim{run}.rds'))
  out
}

psr_sim_fn <- function(n,
                        p,
                        delta = 0.5,
                        compliance_p = 0.75,
                        fn_list, 
                        psfm,
                        run,
                        tmpdir) {
  # browser()
  library(purrr)
  library(tidyr)
  library(dplyr)
  library(clustermq)
  library(synthate)
  library(glue)
  
  bootfn <- function(x, i) {
    # browser()
    new_data <- x[i,]
    
    ps_model <- glm(psfm, family = binomial, data = new_data %>% filter(z == 1))
    new_data <- new_data %>%
      mutate(pr_score = predict(ps_model, newdata = new_data, type = 'response'),
             ps_grp = Hmisc::cut2(pr_score, g = 5),
             pi_c = mean(pr_score)) %>%
      group_by(ps_grp) %>%
      mutate(pi_cj = mean(pr_score)) %>%
      ungroup
    
    
    out <- estimate_ates(new_data,  fn_list)
    as.matrix(out)
  }
  
  set.seed(run)
  print(glue('Simulation {run} for n = {n}, p = {p}, compliance_p = {compliance_p}, delta = {delta}.'))
  sim_data <- 
    generate_data_psr(n, p, delta, compliance_p)
  train_data <- sim_data %>% sample_frac(0.5)
  test_data <- sim_data %>% anti_join(train_data)
  
  #------------------------------
  ## first, all analysis on full data
  #-----------------------------
  ps_model <- glm(psfm, family = binomial, data = sim_data %>% filter(z == 1))
  
  sim_data <- sim_data %>%
    mutate(pr_score = predict(ps_model, newdata = sim_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  full_caces <- estimate_ates(sim_data, fn_list)
  b_theta <- boot::boot(sim_data, bootfn, R = 200)$t
  boot_theta_s <- combine_estimators(ests = full_caces, boot_ests = b_theta)
  full_synth_ests <- boot_theta_s$ate_res %>%
    group_by(theta_0, shrunk, synthetic) %>%
    transmute(estimate = ate,
              var = var,
              sample = 'full') %>%
    ungroup
  
  #----------------------------
  ## now use train data to fit model and combine on test data
  #---------------------------
  ps_model_1 <- glm(psfm, family = binomial, data = train_data %>% filter(z == 1))
  train_data <- train_data %>%
    mutate(pr_score = predict(ps_model_1, newdata = train_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  test_data <- test_data %>%
    mutate(pr_score = predict(ps_model_1, newdata = test_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  train_caces_1 <- estimate_ates(train_data, fn_list)
  test_caces_1 <- estimate_ates(test_data, fn_list)
  
  # b_theta_1 <- boot::boot(train_data, bootfn, R = 200)$t
  boot_theta_s_1 <- combine_estimators(ests = train_caces_1, boot_ests = b_theta)
  
  
  # bhat_1 <- boot_theta_s_1$b_res %>% filter(!shrunk) %>% select(b) %>% unlist
  # holdout_est_1 <- unlist(test_caces_1) %*% bhat_1
  
  test_caces_ds <- test_caces_1 %>%
    gather(est, cace)
  holdout_ests_1 <- test_caces_ds %>%
    inner_join(boot_theta_s_1$b_res) %>%
    group_by(theta_0, shrunk) %>%
    summarise(estimate = sum(cace*b))%>%
    mutate(sample = 'CV',
           synthetic = TRUE)
  
  #-------------------------------
  ## flip the script
  #-----------------------------
  ps_model_2 <- glm(psfm, family = binomial, data = test_data %>% filter(z == 1))
  train_data <- train_data %>%
    mutate(pr_score = predict(ps_model_2, newdata = train_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  test_data <- test_data %>%
    mutate(pr_score = predict(ps_model_2, newdata = test_data, type = 'response'),
           ps_grp = Hmisc::cut2(pr_score, g = 5),
           pi_c = mean(pr_score)) %>%
    group_by(ps_grp) %>%
    mutate(pi_cj = mean(pr_score)) %>%
    ungroup
  train_caces_2 <- estimate_ates(train_data, fn_list)
  test_caces_2 <- estimate_ates(test_data, fn_list)
  
  # b_theta_2 <- boot::boot(test_data, bootfn, R = 200)$t
  boot_theta_s_2 <- combine_estimators(ests = test_caces_2, boot_ests = b_theta)
  
  # bhat_2 <- boot_theta_s_2$b_res %>% filter(!shrunk) %>% select(b) %>% unlist
  # holdout_est_2 <- unlist(train_caces_2) %*% bhat_2
  train_caces_ds <- train_caces_2 %>%
    gather(est, cace)
  holdout_ests_2 <- train_caces_ds %>%
    inner_join(boot_theta_s_2$b_res) %>%
    group_by(theta_0, shrunk) %>%
    summarise(estimate = sum(cace*b)) %>%
    mutate(sample = 'CV',
           synthetic = TRUE) %>%
    ungroup
  
  holdout_ests <- holdout_ests_1 %>%
    full_join(holdout_ests_2) %>%
    group_by(theta_0, shrunk, sample, synthetic) %>%
    summarise(estimate = mean(estimate)) %>%
    ungroup
  
  long_caces  <- full_caces %>%
    gather(theta_0, estimate) %>%
    mutate(var = diag(cov(b_theta)),
           shrunk = FALSE,
           sample = 'full',
           synthetic = FALSE) %>%
    full_join(full_synth_ests) %>%
    full_join(holdout_ests)
  
  out <- long_caces %>%
    mutate(n = n,
           p = p,
           compliance_p = compliance_p,
           delta = delta)
  saveRDS(out, 
          glue('{tmpdir}psrg-compliance-res-n{n}-p{p}-delta{delta}-compliance_p{compliance_p}-sim{run}.rds'))
  out
}
