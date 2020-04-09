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
#' @importFrom mgcv gam
#' @importFrom tibble tibble
#' 

cace <- function(ds, specific_cace, ...) {
  specific_cace(ds, ...)
}

#' @rdname cace
#' @export 
iv_fn <- function(ds) {
  # browser()
  iv_est <- (ds %>%
               group_by(z) %>%
               summarise(ybar = mean(y/pi_c)) %>%
               ungroup %>%
               select(ybar) %>%
               unlist %>%
               diff
  )
  data.frame(iv_est = iv_est)
}

#' @rdname cace
#' @export 
ivw_fn <- function(ds) {
  # browser()
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
tsls_fn <- function(ds, sfm = NULL, yfm = NULL) {
  if (is.null(sfm)) {
    tsls_m1 <- lm(s ~ x + z, data = ds)
  } else {
    tsls_m1 <- lm(sfm, data = ds)
  }
  ds <- ds %>%
    mutate(shat = predict(tsls_m1, newdata = ds))
  if (is.null(yfm)) {
    fit <- lm(y ~ shat + x, data = ds)
  } else {
    fit <- lm(yfm, data = ds)
  }
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
atregr_fn <- function(ds, fm = NULL) {
  # browser()
  # fit <- mgcv::gam(y ~ s + s(pr_score), data = ds)
  if (is.null(fm)) {
    fit <- lm(y ~ s + x, data = ds)
  } else {
    fit <- lm(fm, data = ds)
  }
  
  atregr_est <- coef(fit)[names(coef(fit)) == 's']
  data.frame(atregr_est)
}
#' @rdname cace
#' @export 
ppregr_fn <- function(ds, fm = fm) {
  if (is.null(fm)) {
    fit <- lm(y ~ s + x, data = filter(ds, z == s))
  } else {
    fit <- lm(fm, data = filter(ds, z == s))
  }
  
  ppregr_est <- coef(fit)[names(coef(fit)) == 's']
  data.frame(ppregr_est)
}
#' @rdname cace
#' @export 
ivs_fn <- function(ds) {
  strat_ests <- ds %>%
    group_by(ps_grp, pi_cj, pi_c) %>%
    do(iv_fn(.)) %>%
    ungroup %>%
    filter(iv_est < Inf,
           iv_est > -Inf,
           !is.na(iv_est),
           !is.nan(iv_est)) %>%
    summarise(ivs_est = mean(iv_est*pi_cj/pi_c))
}
#' @rdname cace
#' @export 
ats_fn <- function(ds, fm = NULL) {
  ds %>%
    group_by(ps_grp, pi_cj, pi_c) %>%
    do(atregr_fn(., fm = fm)) %>%
    ungroup %>%
    filter(atregr_est < Inf,
           atregr_est > -Inf,
           !is.na(atregr_est),
           !is.nan(atregr_est)) %>%
    summarise(ats_est = mean(atregr_est*pi_cj/pi_c))
}
#' @rdname cace
#' @export 
pps_fn <- function(ds, fm = NULL) {
  ds %>%
    group_by(ps_grp, pi_cj, pi_c) %>%
    do(ppregr_fn(., fm = fm)) %>%
    ungroup %>%
    filter(ppregr_est < Inf,
           ppregr_est > -Inf,
           !is.na(ppregr_est),
           !is.nan(ppregr_est)) %>%
    summarise(pps_est = mean(ppregr_est*pi_cj/pi_c))
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
ipw_fn <- function(ds, p = NULL) {
  if (is.null(p)) {
    ipw_est <- with(ds, 
                    PSPS_SM_weighting(Z = z,
                                      D = s,
                                      X = x,
                                      Y = y))$CACE
  } else {
    ipw_est <- with(ds, 
                    PSPS_SM_weighting(Z = z,
                                      D = s,
                                      X = as.matrix(dplyr::select(ds, starts_with('x'))),
                                      Y = y))$CACE
  }
  
  tibble::tibble(ipw_est = ipw_est)
}
#' @rdname cace
#' @export 
ipw_regr_fn <- function(ds, p = NULL) {
  if (is.null(p)) {
    ipw_regr_est <- with(ds, 
                         PSPS_SM_weighting(Z = z,
                                           D = s,
                                           X = x,
                                           Y = y))$CACE.reg
  } else {
    ipw_regr_est <- with(ds, 
                         PSPS_SM_weighting(Z = z,
                                           D = s,
                                           X = as.matrix(dplyr::select(ds, starts_with('x'))),
                                           Y = y))$CACE.reg
  }
  tibble::tibble(ipw_regr_est = ipw_regr_est)
}
#' @rdname cace
#' @export 
ipws_fn <- function(ds, p = NULL) {
  ds %>%
    group_by(ps_grp, pi_cj, pi_c) %>%
    do(ipw_fn(., p = p)) %>%
    ungroup %>%
    filter(ipw_est < Inf,
           ipw_est > -Inf,
           !is.na(ipw_est),
           !is.nan(ipw_est)) %>%
    summarise(ipws_est = mean(ipw_est*pi_cj/pi_c))
}
#' @rdname cace
#' @export 
ipwrs_fn <- function(ds, p = NULL) {
  ds %>%
    group_by(ps_grp, pi_cj, pi_c) %>%
    do(ipw_regr_fn(., p = p)) %>%
    ungroup %>%
    filter(ipw_regr_est < Inf,
           ipw_regr_est > -Inf,
           !is.na(ipw_regr_est),
           !is.nan(ipw_regr_est)) %>%
    summarise(ipwrs_est = mean(ipw_regr_est*pi_cj/pi_c))
}
PSPS_SM_weighting = function(Z, D, X, Y, ep = 1) {
  # browser()
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
  if (length(index11) < ncol(X)) return(
    list(CACE = NA, #NACE = NACE, 
         CACE.reg = NA)
  )
  
  #weights for regression estimator
  wc = ps.score.c[Z==0]/pr.compliers
  wn = ps.score.n[Z==0]/pr.nevertakers
  
  #model assisted regression estimator 
  r1c = lm(Y[index11] ~ 0 + X[index11, ])$coef
  # r1n = lm(Y[index10] ~ 0 + X[index10, ])$coef
  r0c = lm(Y[Z==0] ~ 0 + X[Z==0, ], weights = wc)$coef
  # r0n = lm(Y[Z==0] ~ 0 + X[Z==0, ], weights = wn)$coef
  
  #weighted outcomes
  weighted.Y.c = Y[Z==0]*wc
  # weighted.Y.n = Y[Z==0]*wn
  
  #CACE and NACE
  CACE = mean(Y[index11]) - mean(weighted.Y.c)
  # NACE = mean(Y[index10]) - mean(weighted.Y.n)
  
  #weighted outcomes for regression estimator
  weighted.Y1c = Y[index11]-X[index11, ]%*%r1c
  weighted.Y0c = (Y[Z==0]-X[Z==0, ]%*%r0c)*wc
  # weighted.Y1n = Y[index10]-X[index10, ]%*%r1n
  # weighted.Y0n = (Y[Z==0]-X[Z==0, ]%*%r0n)*wn
  weighted.rc = rbind(X[index11, ], X[Z==0, ]*wc) %*% (r1c - r0c)
  # weighted.rn = rbind(X[index10, ], X[Z==0, ]*wn) %*% (r1n - r0n)
  
  #CACE and NACE regression estimates
  CACE.reg = mean(weighted.Y1c) - mean(weighted.Y0c) + mean(weighted.rc)
  # NACE.reg = mean(weighted.Y1n) - mean(weighted.Y0n) + mean(weighted.rn)
  
  #results
  ACE = list(CACE = CACE, #NACE = NACE, 
             CACE.reg = CACE.reg)#, NACE.reg = NACE.reg)
  return(ACE)
}

#' @rdname cace
#' @export 
get_weighted_y_c = function(Z, D, X, Y) {
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
  
  #the probability of compliers and never-takers
  pr.compliers   = sum(Z*D)/sum(Z)
  
  #weights for regression estimator
  wc = ps.score.c[Z==0]/pr.compliers
  
  #weighted outcomes
  weighted.Y.c = Y[Z==0]*wc
  weighted.Y.c
}
