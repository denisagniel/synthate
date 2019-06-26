#' Functions to estimate the average treatment effect
#'
#' @param this_data a data frame with value \code{y} for the outcome, \code{d} for the treatment, and \code{w} for a case weight.
#' @param ... additional arguments to be passed indicating, for example, the outcome formula, propensity score formula, or other things on a case by case basis.
#'
#' @return a 1 x 1 data frame containing the named ATE
#' @export 
#' @importFrom randomForest randomForest
#' @importFrom stringr str_c
#' @importFrom rpart rpart rpart.control
#' @importFrom Matching Match
#' @importFrom dplyr select select_ one_of
#' @importFrom grf causal_forest average_treatment_effect
#' @importFrom ATE ATE
#' @importFrom magrittr "%>%"
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
#' this_data <- estimate_scores(this_data, outcome_fm = outcome_fm,
#'                              ps_fm = ps_fm,
#'                              ps_fam = ps_fam,
#'                              outcome_fam = outcome_fam)
#' thetahat <- with(gen_mod,
#'                  estimate_ates(this_data,
#'                           ate_list,
#'                           cov_ids = cov_ids,
#'                           outcome_fm = stringr::str_c('d + ', outcome_fm),
#'                           outcome_fam = outcome_fam))
#' thetahat
#' 

ate <- function(this_data, specific_ate, ...) {
  specific_ate(this_data, ...)
}

#' @rdname ate
#' @export 
naive_ate <- function(this_data, ...) {
  this_data %>%
    summarise(ate_naive = sum(y*d*w)/sum(d*w) - sum(y*(1-d)*w)/sum((1-d)*w))
}

#' @rdname ate
#' @export
ipw1_ate <- function(this_data, ...) {
  this_data %>%
    summarise(ate_ipw_1 = mean(y/pi*d*w) - mean(y/pi*(1-d)*w))
}

#' @rdname ate
#' @export
ipw2_ate <- function(this_data, ...) {
  this_data %>%
    summarise(ate_ipw_2 = sum(y/pi*d*w)/sum(d/pi*w) - sum(y/pi*(1-d)*w)/sum((1-d)/pi*w))
}

#' @rdname ate
#' @export
ipw3_ate <- function(this_data, ...) {
  this_data %>%
    mutate(c1 = sum((d - ps)/ps)/sum((d-ps)^2/ps^2),
           c0 = -sum((d - ps)/(1-ps))/sum((d-ps)^2/(1-ps)^2),
           a = d/ps*(1-c1/ps),
           b = (1-d)/(1-ps)*(1-c0/(1-ps))) %>%
    summarise(ate_ipw_3 = sum(y*a*w)/sum(a*w) - sum(y*b*w)/sum(b*w))
}

#' @rdname ate
#' @export
strat_ate <- function(this_data, ...) {
  strata_summary <- this_data %>% group_by(d, ps_strat) %>%
    summarise(mu = sum(y*w)/sum(w))
  strata_n <- this_data %>% group_by(ps_strat) %>%
    summarise(n_s = sum(w))
  strata_summary %>% inner_join(strata_n, by = 'ps_strat') %>%
    group_by(d) %>%
    summarise(mu = sum(mu*n_s)/sum(n_s)) %>%
    ungroup %>%
    summarise(ate_strat = diff(mu))
}

#' @rdname ate
#' @export
rf_ate <- function(this_data, outcome_fm, ...) {
  # browser()
  outcome_fm <- as.formula(stringr::str_c('y ~ d + ', outcome_fm))
  fit_0 <- randomForest::randomForest(outcome_fm, data = this_data %>% filter(d == 0))
  fit_1 <- randomForest::randomForest(outcome_fm, data = this_data %>% filter(d == 1))
  this_data <- this_data %>% mutate(
    yhat0 = predict(fit_0, newdata = this_data %>% mutate(d = 0)),
    yhat1 = predict(fit_1, newdata = this_data %>% mutate(d = 1))
  )

  this_data %>% summarise(rf_ate =  sum(w*(yhat1 - yhat0))/sum(w))
}

#' @rdname ate
#' @export
cart_ate <- function(this_data, outcome_fm, ...) {
  # browser()
  outcome_fm <- as.formula(stringr::str_c('y ~ d + ', str_replace(outcome_fm, '\\*', '\\+')))
  fit_0 <- rpart::rpart(outcome_fm, data = this_data %>% filter(d == 0),
                        # weights = 1/pi,
                         control = rpart::rpart.control(
                           minbucket = 1,
                           cp = 1e-30,
                           maxdepth = 30,
                           xval = 1
                         )
  )
  fit_1 <- rpart::rpart(outcome_fm, data = this_data %>% filter(d == 1),
                        # weights = 1/pi,
                        control = rpart::rpart.control(
                          minbucket = 1,
                          cp = 1e-30,
                          maxdepth = 30,
                          xval = 1
                        )
  )
  # fit <- rpart::rpart(outcome_fm, data = this_data,
  #                                           control = rpart::rpart.control(
  #                                             minbucket = 5,
  #                                             cp = 1e-30,
  #                                             maxdepth = 30,
  #                                             xval = 1
  #                                           ))
  this_data <- this_data %>% mutate(
    yhat0 = predict(fit_0, newdata = this_data %>% mutate(d = 0)),
    yhat1 = predict(fit_1, newdata = this_data %>% mutate(d = 1)),
    # mu_0 = -2.8 + 1.5*(x1 > 1) + 2*(x2 < -1) + 3*(x3 > 1)*(x1 < 1) + 4*(x4 < 1) + 5*(x5 > 0) +
      # 5*v1 + 4*v2 + 3*v3 + 2*v4 + 1.5*v5,
    # mu_1 = mu_0 + 1
  )

  this_data %>% summarise(cart_ate =  sum(w*(yhat1 - yhat0))/sum(w))
}

#' @rdname ate
#' @export
pscart_ate <- function(this_data, ps_fm, ...) {
  browser()
  this_fm <- as.formula(stringr::str_c('d ~ ', ps_fm))
  fit <- rpart::rpart(this_fm, data = this_data,
                        control = rpart::rpart.control(
                          minbucket = 5,
                          cp = 1e-30,
                          maxdepth = 30,
                          xval = 1
                        ))
  party_fit <- as.party(fit)
  fit_node <- fitted(party_fit)
  this_data <- this_data %>% mutate(pscart_node = fit_node[,1])

  this_data %>% group_by(pscart_node) %>%
    summarise(ydiff =  case_when(
      all(d == 1) | all(d == 0) ~ as.double(NA),
      TRUE ~ sum((y*w)[d == 1] - (y*w)[d==0])/sum(w))
    ) %>%
    ungroup %>% summarise(pscart_ate = mean(ydiff, na.rm = TRUE))
}

#' @rdname ate
#' @export
strat_regr_ate <- function(this_data, outcome_fm, outcome_fam = gaussian, ...) {
  # browser()
  outcome_fm <- as.formula(stringr::str_c('y ~ ', outcome_fm))
  stratified_diffs <- lapply(unique(this_data$ps_strat), function(i) {
    # browse r()
    this_sub <- this_data %>% filter(ps_strat == i)
    ini_fit <- lm(outcome_fm, data = this_sub, weights = w)
    this_sub <- this_sub %>% mutate(yh_tmp = predict(ini_fit))
    if (!identical(outcome_fam, gaussian)) {
      fit <- glm(outcome_fm, start = ini_fit$coefficients,
                 etastart = yh_tmp,
                 data = this_sub,
                 family = outcome_fam)
      # dsgn <- survey::svydesign(ids = ~1,
      #                           weights = ~w, data = this_data,
      #                           formula = outcome_fm)
      # fit <- Dmisc::myTry(survey::svyglm(outcome_fm,
      #                       data = this_data,
      #                       family = outcome_fam,
      #                       start = ini_fit$coefficients,
      #                       design = dsgn))
      if (class(fit) == 'try-error') return(data.frame(te = NA, ps_strat = i))
    } else fit <- ini_fit
    y0 <- predict(fit, newdata = this_sub %>% mutate(d = 0), type = 'response')
    y1 <- predict(fit, newdata = this_sub %>% mutate(d = 1), type = 'response')
    this_sub %>% summarise(te = sum(w*(y1 - y0))/sum(w), ps_strat = i)
  }) %>% bind_rows
  # browser()
  strata_n <- this_data %>% group_by(ps_strat) %>%
    summarise(n_s = n())
  stratified_diffs %>% inner_join(strata_n, by = 'ps_strat') %>%
    summarise(ate_sr = sum(te*n_s, na.rm = TRUE)/sum(n_s))
}

generic_match_ate <- function(this_data, matching_vars,
                          M = 5, nm, ...) {
  # browser()
  match <- with(this_data, Matching::Match(Y = y, Tr = d, X = matching_vars,
                                           estimand = 'ATE', version = 'standard',
                                           ties = FALSE, M = M, ...))
  if (!is.na(match)) {
    df <- data.frame(match$est)
  } else df <- data.frame(NA)
  colnames(df) <- nm
  df
}

#' @rdname ate
#' @export
match_ps_ate <- function(this_data, ...) {
  if (is.null(this_data$ps_K)) {
    generic_match_ate(this_data, this_data$ps, nm = 'match_ps')
  } else {
    this_data %>% summarise(match_ps_ate = sum(y*w*ps_K*d)/sum(w*ps_K*d) -
                sum(y*w*ps_K*(1-d))/sum(w*ps_K*(1-d)))
  }

}
#' @rdname ate
#' @export
match_prog_ate <- function(this_data, ...) {
  if (is.null(this_data$prog_K)) {
    generic_match_ate(this_data, this_data$prog_score, nm = 'match_prog')
  } else {
    this_data %>% summarise(match_prog_ate = sum(y*w*prog_K*d)/sum(w*prog_K*d) -
                sum(y*w*prog_K*(1-d))/sum(w*prog_K*(1-d)))
  }
}
#' @rdname ate
#' @export
match_both_ate <- function(this_data, ...) {
  if (is.null(this_data$both_K)) {
    generic_match_ate(this_data, cbind(this_data$ps, this_data$prog_score), nm = 'match_both')
  } else {
    this_data %>% summarise(match_both_ate = sum(y*w*both_K*d)/sum(w*both_K*d) -
                              sum(y*w*both_K*(1-d))/sum(w*both_K*(1-d)))
  }
}

#' @rdname ate
#' @export
caliper_ps_ate <- function(this_data, ...) {
  if (is.null(this_data$cal_ps_K)) {
    generic_match_ate(this_data, this_data$ps, nm = 'cal_match_ps',
                      caliper = 0.05)
  } else {
    this_data %>%
      summarise(
        cal_match_ps = sum(y*w*cal_ps_K*d, na.rm = TRUE) /
          sum(w*cal_ps_K*d, na.rm = TRUE) -
          sum(y*w*cal_ps_K*(1-d), na.rm = TRUE) /
          sum(w*cal_ps_K*(1-d), na.rm = TRUE))
  }

}
#' @rdname ate
#' @export
caliper_prog_ate <- function(this_data, ...) {
  if (is.null(this_data$cal_prog_K)) {
    generic_match_ate(this_data, this_data$prog_score, nm = 'cal_match_prog',
                      caliper = 0.05)
  } else {
    this_data %>%
      summarise(
        cal_match_prog = sum(y*w*cal_prog_K*d, na.rm = TRUE) /
          sum(w*cal_prog_K*d, na.rm = TRUE) -
          sum(y*w*cal_prog_K*(1-d), na.rm = TRUE) /
          sum(w*cal_prog_K*(1-d), na.rm = TRUE))
  }

}
#' @rdname ate
#' @export
caliper_both_ate <- function(this_data, ...) {
  if (is.null(this_data$cal_both_K)) {
    generic_match_ate(this_data, cbind(this_data$ps, this_data$prog_score),
                      nm = 'cal_match_both', caliper = 0.2)
    } else {
    this_data %>%
        summarise(
          cal_match_both = sum(y*w*cal_both_K*d, na.rm = TRUE) /
            sum(w*cal_both_K*d, na.rm = TRUE) -
            sum(y*w*cal_both_K*(1-d), na.rm = TRUE) /
            sum(w*cal_both_K*(1-d), na.rm = TRUE))
  }

}
#' @rdname ate
#' @export
caliper_nn_ate <- function(this_data, ...) {
  if (is.null(this_data$cal_nn_K)) {
    generic_match_ate(this_data, this_data$ps, nm = 'cal_match_nn',
                      caliper = 0.05, replace = FALSE, M = 1)
  } else {
    this_data %>%
      summarise(cal_match_nn = sum(y*w*cal_nn_K*d, na.rm = TRUE) /
                  sum(w*cal_nn_K*d, na.rm = TRUE) -
                  sum(y*w*cal_nn_K*(1-d), na.rm = TRUE) /
                  sum(w*cal_nn_K*(1-d), na.rm = TRUE))
  }

}

#' @rdname ate
#' @export
caliper_logit_ate <- function(this_data, ...) {
  if (is.null(this_data$cal_logit_K)) {
    generic_match_ate(this_data, log(this_data$ps/(1-this_data$ps)), nm = 'cal_match_logit',
                      caliper = 0.05, replace = FALSE, M = 1)
  } else {
    this_data %>%
      summarise(cal_match_logit =
                  sum(y*w*cal_logit_K*d, na.rm = TRUE) /
                  sum(w*cal_logit_K*d, na.rm = TRUE) -
                  sum(y*w*cal_logit_K*(1-d), na.rm = TRUE) /
                  sum(w*cal_logit_K*(1-d), na.rm = TRUE))
  }

}

#' @rdname ate
#' @export
regr_ate <- function(this_data, outcome_fm, outcome_fam = gaussian, ...) {
  # browser()
  outcome_fm <- as.formula(stringr::str_c('y ~ ', outcome_fm))
  lm_ini <- lm(outcome_fm, data = this_data, weights = w)
  b_ini <- coef(lm_ini)
  yh_ini <- predict(lm_ini)
  if (!identical(outcome_fam, gaussian)) {
    outcome_fit <- glm(outcome_fm, data = this_data, family = outcome_fam,
                  start = coef(lm_ini),
                  etastart = abs(predict(lm_ini)))
    # if (Dmisc::isErr(outcome_fit)) browser()
    # dsgn <- survey::svydesign(ids = ~1,
    #                           weights = ~w, data = this_data,
    #                           formula = outcome_fm)
    # # lmic <- lm_ini$coefficients
    # outcome_fit <- survey::svyglm(outcome_fm,
    #                       data = this_data,
    #                       family = outcome_fam,
    #                       # start = lmic,
    #                       design = dsgn,
    #                       etastart = yh_ini)
  } else outcome_fit <- lm_ini
  this_data <- this_data %>%
    mutate(yhat_1 = predict(outcome_fit, newdata = this_data %>% mutate(d = 1),
                            type = 'response'),
           yhat_0 = predict(outcome_fit, newdata = this_data %>% mutate(d = 0),
                            type = 'response'))
  this_data %>% summarise(ate_regr = sum(w*(yhat_1 - yhat_0))/sum(w))
}

#' @rdname ate
#' @export
dr_ate <- function(this_data, outcome_fm, outcome_fam = gaussian, ...) {
  # browser()
  outcome_fm <- as.formula(stringr::str_c('y ~ ', outcome_fm))
  lm_ini <- lm(outcome_fm, data = this_data, weights = w)
  if (!identical(outcome_fam, gaussian)) {
    outcome_fit <- glm(outcome_fm, data = this_data, family = outcome_fam,
                       start = coef(lm_ini), etastart = abs(predict(lm_ini)))
    # if (Dmisc::isErr(outcome_fit)) browser()
    # dsgn <- survey::svydesign(ids = ~1,
    #                           weights = ~w, data = this_data,
    #                           formula = outcome_fm)
    # outcome_fit <- survey::svyglm(outcome_fm,
    #                       data = this_data,
    #                       family = outcome_fam,
    #                       design = dsgn,
    #                        etastart = predict(lm_ini))
  } else outcome_fit <- lm_ini
  this_data <- this_data %>%
    mutate(yhat_1 = predict(outcome_fit, newdata = this_data %>% mutate(d = 1),
                            type = 'response'),
           yhat_0 = predict(outcome_fit, newdata = this_data %>% mutate(d = 0),
                            type = 'response'))
  this_data %>%
    summarise(ate_dr = sum(w*(d*y/ps - (d - ps)/ps*yhat_1))/sum(w) -
                sum(w*((1-d)*y/(1-ps) + (d - ps)/(1-ps)*yhat_0))/sum(w))
}

#' @rdname ate
#' @export
bal_ate <- function(this_data, cov_ids, ...) {
  # browser()
  this_covm <- this_data %>% dplyr::select(dplyr::one_of(cov_ids)) %>% as.matrix
  ate_bal <- try(ATE::ATE(Y = this_data$y, Ti = this_data$d,
                          X = this_covm), silent = TRUE)
  if (class(ate_bal) == 'try-error') {
    data.frame(ate_bal = NA)
  } else data.frame(ate_bal = ate_bal$est[3])
}

#' @rdname ate
#' @export
highdim_bal_ate <- function(this_data, cov_ids, ...) {
  # browser()
  this_covm <- this_data %>% dplyr::select(dplyr::one_of(cov_ids)) %>% as.matrix
  bhd_fit <- try(balanceHD::residualBalance.ate(X = this_covm,
                                                Y = this_data$y,
                                                W = this_data$d), silent = TRUE)
  if (class(bhd_fit) == 'try-error') {
    data.frame(ate_bhd = NA)
  } else data.frame(ate_bhd = bhd_fit)
}

#' @rdname ate
#' @export
cem_ate <- function(this_data, cov_ids, ...) {
  # browser()
  bal_data <- this_data %>% dplyr::select(d, one_of(cov_ids))
  bal_1 <- cem::cem(treatment = 'd', data = bal_data)
  if (class(bal_1) == 'try-error') {
    data.frame(ate_cem = NA)
  } else {
    att_1 <- cem::att(bal_1, y ~ d, data = this_data)
    bal_0 <- try(cem::cem('d', data = bal_data, baseline.group = '0'), silent = TRUE)
    if (class(bal_0) == 'try-error') {
      data.frame(ate_cem = NA)
    } else att_0 <- cem::att(bal_0, y ~ d, data = this_data)
    data.frame(ate_cem = (att_1$att.model[1,2]*sum(this_data$d) +
                            att_0$att.model[1,2]*sum(1-this_data$d))/nrow(this_data))
  }
}

#' @rdname ate
#' @export
grf_ate <- function(this_data, cov_ids, ...) {
  x <- this_data %>% dplyr::select_(.dots = cov_ids) %>% as.matrix
  grf_fit <- grf::causal_forest(X = x, W = this_data$d, Y = this_data$y)
  data.frame(grf_ate = grf::average_treatment_effect(grf_fit)[1])
}


het_trt_grf <- function(this_data, cov_ids, ...) {
  x <- this_data %>% dplyr::select_(.dots = cov_ids) %>% as.matrix
  grf_fit <- grf::causal_forest(X = x, W = this_data$d, Y = this_data$y)
  predict(grf_fit)$predictions
}
