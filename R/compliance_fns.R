#' Functions to estimate the complier average causal effect or local average treatment effect
#'
#' @param ds a data frame with value \code{y} for the outcome, \code{z} for randomization/instrument, \code{s} for treatment receipt, and covariate \code{x}.
#' @param ... additional arguments to be passed.
#'
#' @return a 1 x 1 data frame containing the named ATE
#' @export 
#'
#' @examples
#' 
#' 
#' 

cace <- function(ds, specific_cace, ...) {
  specific_cace(ds, ...)
}

#' @rdname cace
#' @export 
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
#' @rdname cace
#' @export 
ivw_fn <- function(ds) {
  # grf_fit <- with(ds, grf::causal_forest(Y = y, X = x, W = z))
  ds_0 <- ds %>% mutate(z = 0)
  ds_1 <- ds %>% mutate(z = 1)
  grf_fit <- mgcv::gam(y ~ s(x, by = z) + s(x), data = ds)
  delta <- predict(grf_fit, newdata = ds_1) - predict(grf_fit, newdata = ds_0)
  data.frame(ivw_fn = sum(delta)/sum(ds$pr_score))
}
#' @rdname cace
#' @export 
at_fn <- function(ds) {
  ybar_s1 <- ds %>% filter(s == 1) %>%
    summarise(mean(y)) %>% unlist
  ybar_s0 <- ds %>% filter(s == 0) %>%
    summarise(mean(y)) %>% unlist

  at_est <- ybar_s1 - ybar_s0
  data.frame(at_est = at_est)
}
#' @rdname cace
#' @export 
pp_fn <- function(ds) {
  ybar_pp1 <- ds %>% filter(s == 1, z == 1) %>%
    summarise(mean(y)) %>% unlist
  ybar_pp0 <- ds %>% filter(s == 0, z == 0) %>%
    summarise(mean(y)) %>% unlist

  pp_est <- ybar_pp1 - ybar_pp0
  data.frame(pp_est = pp_est)
}
#' @rdname cace
#' @export 
tsls_fn <- function(ds) {
  tsls_m1 <- lm(s ~ x + z, data = ds)
  ds <- ds %>%
    mutate(shat = predict(tsls_m1, newdata = ds))
  fit <- lm(y ~ shat + x, data = ds)
  tsls_est <- coef(fit)[2]
  data.frame(tsls_est = tsls_est)
}
#' @rdname cace
#' @export 
itt_regr_fn <- function(ds) {
  # browser()
  # fit <- mgcv::gam(y ~ z + s(pr_score), data = ds)
  fit <- lm(y ~ z + x, data = ds)
  regr_est <- coef(fit)[2]
  data.frame(regr_est)
}
#' @rdname cace
#' @export 
atregr_fn <- function(ds) {
  # browser()
  # fit <- mgcv::gam(y ~ s + s(pr_score), data = ds)
  fit <- lm(y ~ s + x, data = ds)
  atregr_est <- coef(fit)[2]
  data.frame(atregr_est)
}
#' @rdname cace
#' @export 
ppregr_fn <- function(ds) {
  # browser()
  # fit <- mgcv::gam(y ~ s + s(pr_score), data = ds %>% filter(z == s))
  fit <- lm(y ~ s + x, data = ds %>% filter(z == s))
  ppregr_est <- coef(fit)[2]
  data.frame(ppregr_est)
}
#' @rdname cace
#' @export 
ivs_fn <- function(ds) {
  ds %>%
    group_by(ps_grp) %>%
    do(iv_fn(.)) %>%
    ungroup %>%
    summarise(ivs_est = mean(iv_est))
}
#' @rdname cace
#' @export 
ats_fn <- function(ds) {
  ds %>%
    group_by(ps_grp) %>%
    do(at_fn(.)) %>%
    ungroup %>%
    summarise(ats_est = mean(at_est))
}
#' @rdname cace
#' @export 
pps_fn <- function(ds) {
  ds %>%
    group_by(ps_grp) %>%
    do(pp_fn(.)) %>%
    ungroup %>%
    summarise(pps_est = mean(pp_est))
}

# ipw_fn <- function(ds) {
#   ds %>%
#     mutate(
#       psz = mean(s[z == 1]),
#       ipw_w = case_when(
#         z == 1 ~ 1,
#         z == 0 ~ -pr_score/psz
#       )
#   ) %>%
#     summarise(ipw_est = mean(y*ipw_w))
# }
#' @rdname cace
#' @export 
ipw_fn <- function(ds) {
  with(ds, 
       PSPS_SM_weighting(Z = z,
                         D = s,
                         X = x,
                         Y = y))$CACE
}
#' @rdname cace
#' @export 
ipw_regr_fn <- function(ds) {
  with(ds, 
       PSPS_SM_weighting(Z = z,
                         D = s,
                         X = x,
                         Y = y))$CACE.reg
}

PSPS_SM_weighting = function(Z, D, X, Y, ep = 1) {
  #augment the design X
  N = length(Z) 
  X = cbind(rep(1, N), X)
  
  #estimate the weights
  D1 = D[Z==1]
  X1 = X[Z==1, ]
  data1 = cbind(D1, X1)
  data1 = as.data.frame(data1)
  
  #logistic regression
  logit.treat = glm(D1 ~ 0 +., data = data1, family = binomial(link = logit))
  beta = coef(logit.treat)  
  
  #the predicted propensity score
  #for compliers
  ps.score.c = 	1/(1 + exp(-X%*%beta))
  #for never-takers
  ps.score.n = 1 - ps.score.c
  
  #adding sensitivity parameter into score model
  ps.score.c = ep*ps.score.c/(ep*ps.score.c + ps.score.n)
  ps.score.n = ps.score.n/(ep*ps.score.c + ps.score.n)
  
  #the probability of compliers and never-takers
  pr.compliers   = sum(Z*D)/sum(Z)
  pr.nevertakers = 1 - pr.compliers
  
  #indices
  index11 = (1:N)[Z==1&D==1]
  index10 = (1:N)[Z==1&D==0]
  index01 = (1:N)[Z==0&D==1]
  index00 = (1:N)[Z==0&D==0]
  
  #weights for regression estimator
  wc = ps.score.c[Z==0]/pr.compliers
  wn = ps.score.n[Z==0]/pr.nevertakers
  
  #model assisted regression estimator 
  r1c = lm(Y[index11] ~ 0 + X[index11, ])$coef
  r1n = lm(Y[index10] ~ 0 + X[index10, ])$coef
  r0c = lm(Y[Z==0] ~ 0 + X[Z==0, ], weights = wc)$coef
  r0n = lm(Y[Z==0] ~ 0 + X[Z==0, ], weights = wn)$coef
  
  #weighted outcomes
  weighted.Y.c = Y[Z==0]*wc
  weighted.Y.n = Y[Z==0]*wn
  
  #CACE and NACE
  CACE = mean(Y[index11]) - mean(weighted.Y.c)
  NACE = mean(Y[index10]) - mean(weighted.Y.n)
  
  #weighted outcomes for regression estimator
  weighted.Y1c = Y[index11]-X[index11, ]%*%r1c
  weighted.Y0c = (Y[Z==0]-X[Z==0, ]%*%r0c)*wc
  weighted.Y1n = Y[index10]-X[index10, ]%*%r1n
  weighted.Y0n = (Y[Z==0]-X[Z==0, ]%*%r0n)*wn
  weighted.rc = rbind(X[index11, ], X[Z==0, ]*wc) %*% (r1c - r0c)
  weighted.rn = rbind(X[index10, ], X[Z==0, ]*wn) %*% (r1n - r0n)
  
  #CACE and NACE regression estimates
  CACE.reg = mean(weighted.Y1c) - mean(weighted.Y0c) + mean(weighted.rc)
  NACE.reg = mean(weighted.Y1n) - mean(weighted.Y0n) + mean(weighted.rn)
  
  #results
  ACE = list(CACE = CACE, NACE = NACE, CACE.reg = CACE.reg, NACE.reg = NACE.reg)
  return(ACE)
}

