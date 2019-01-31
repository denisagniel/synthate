
iv_fn <- function(ds) {
  # browser()
  prob_comply <- ds %>% filter(z == 1) %>%
    summarise(mean(c)) %>% unlist
  iv_est <- (ds %>%
      group_by(z) %>%
      summarise(ybar = mean(y)) %>%
      ungroup %>%
      select(ybar) %>%
      unlist %>%
      diff
  )/prob_comply
  data.frame(iv_est = iv_est)
}
#
ivw_fn <- function(ds) {
  # grf_fit <- with(ds, grf::causal_forest(Y = y, X = x, W = z))
  ds_0 <- ds %>% mutate(z = 0)
  ds_1 <- ds %>% mutate(z = 1)
  grf_fit <- mgcv::gam(y ~ s(x, by = z) + s(x), data = ds)
  delta <- predict(grf_fit, newdata = ds_1) - predict(grf_fit, newdata = ds_0)
  data.frame(ivw_fn = sum(delta)/sum(ds$pr_score))
}

at_fn <- function(ds) {
  ybar_s1 <- ds %>% filter(s == 1) %>%
    summarise(mean(y)) %>% unlist
  ybar_s0 <- ds %>% filter(s == 0) %>%
    summarise(mean(y)) %>% unlist

  at_est <- ybar_s1 - ybar_s0
  data.frame(at_est = at_est)
}

pp_fn <- function(ds) {
  ybar_pp1 <- ds %>% filter(s == 1, z == 1) %>%
    summarise(mean(y)) %>% unlist
  ybar_pp0 <- ds %>% filter(s == 0, z == 0) %>%
    summarise(mean(y)) %>% unlist

  pp_est <- ybar_pp1 - ybar_pp0
  data.frame(pp_est = pp_est)
}

tsls_fn <- function(ds) {
  tsls_m1 <- lm(s ~ x + z, data = ds)
  ds <- ds %>%
    mutate(shat = predict(tsls_m1, newdata = ds))
  fit <- lm(y ~ shat + x, data = ds)
  tsls_est <- coef(fit)[2]
  data.frame(tsls_est = tsls_est)
}

itt_regr_fn <- function(ds) {
  # browser()
  # fit <- mgcv::gam(y ~ z + s(pr_score), data = ds)
  fit <- lm(y ~ z + x, data = ds)
  regr_est <- coef(fit)[2]
  data.frame(regr_est)
}

atregr_fn <- function(ds) {
  # browser()
  # fit <- mgcv::gam(y ~ s + s(pr_score), data = ds)
  fit <- lm(y ~ s + x, data = ds)
  atregr_est <- coef(fit)[2]
  data.frame(atregr_est)
}

ppregr_fn <- function(ds) {
  # browser()
  # fit <- mgcv::gam(y ~ s + s(pr_score), data = ds %>% filter(z == s))
  fit <- lm(y ~ s + x, data = ds %>% filter(z == s))
  ppregr_est <- coef(fit)[2]
  data.frame(ppregr_est)
}

ivs_fn <- function(ds) {
  ds %>%
    group_by(ps_grp) %>%
    do(iv_fn(.)) %>%
    ungroup %>%
    summarise(ivs_est = mean(iv_est))
}

ats_fn <- function(ds) {
  ds %>%
    group_by(ps_grp) %>%
    do(at_fn(.)) %>%
    ungroup %>%
    summarise(ats_est = mean(at_est))
}

pps_fn <- function(ds) {
  ds %>%
    group_by(ps_grp) %>%
    do(pp_fn(.)) %>%
    ungroup %>%
    summarise(pps_est = mean(pp_est))
}

ipw_fn <- function(ds) {
  ds %>%
    mutate(
      psz = mean(s[z == 1]),
      ipw_w = case_when(
        z == 1 ~ 1,
        z == 0 ~ -pr_score/psz
      )
  ) %>%
    summarise(ipw_est = mean(y*ipw_w))
}
