naive_omega <- function(this_data, ...) {
  this_data %>%
    transmute(naive_omega = (d/sum(d) - (1-d)/sum((1-d)))/n)
}

ipw1_omega <- function(this_data, ...) {
  this_data %>%
    transmute(ipw1_omega = d/pi/sum(d) - (1-d)/pi/sum(1-d))
}

ipw2_omega <- function(this_data, ...) {
  this_data %>%
    transmute(ipw2_omega = d/pi/sum(d/pi) - (1-d)/pi/sum((1-d)/pi))
}

ipw3_omega <- function(this_data, ...) {
  this_data %>%
    mutate(c1 = sum((d - ps)/ps)/sum((d-ps)^2/ps^2),
       c0 = -sum((d - ps)/(1-ps))/sum((d-ps)^2/(1-ps)^2),
       a = d/ps*(1-c1/ps),
       b = (1-d)/(1-ps)*(1-c0/(1-ps))) %>%
  transmute(ipw3_omega = a/sum(a) - b/sum(b))
}

regr_omega <- function(this_data, outcome_fm, ...) {
  ##--------------------
  ## this assumes OLS
  ##-------------------
  X <- model.matrix(as.formula(stringr::str_c('~ d +', outcome_fm)),
                    data = this_data)
  h <- solve(t(X) %*% X) %*% t(X)
  data.frame(regr_omega = h[2,])
}

# sr_omega <- function(this_data, outcome_fm) {
#   this_data %>% group_by(ps_strat) %>%
#     do(regr_omega(., outcome_fm))
# }

dr_omega <- function(this_data, outcome_fm, ...) {
  X <- model.matrix(as.formula(stringr::str_c('~ d +', outcome_fm)),
                    data = this_data)
  h <- solve(t(X) %*% X) %*% t(X)
  X_1 <- X
  X_1[,2] <- 1
  X_0 <- X
  X_0[,2] <- 0

  h_0 <- X_0 %*% h
  h_1 <- X_1 %*% h
  n <- nrow(this_data)
  this_data %>%
    transmute(dr_omega = n^(-1)*(d/pi - (1-d)/pi - t(h_1) %*% ((d - ps)/ps) -
                t(h_0) %*% ((d - ps)/(1-ps))))
}

match_omega <- function(this_data, matching_vars, nm, M = 5, ...) {
  match_out <-  with(this_data,
                             Matching::Match(Y = y, Tr = d, X = matching_vars,
                                             estimand = 'ATE', version = 'fast',
                                             ties = FALSE, M = M, ...))
  omega <- sapply(1:n, function(i) {
    sum(match_out$index.treated == i)/length(match_out$index.treated) -
      sum(match_out$index.control == i)/length(match_out$index.control)
  })
  tmp <- data.frame(omega)
  colnames(tmp) <- nm
  tmp
}

bal_omega <- function(this_data, covar_ids, ...) {
  this_covm <- this_data %>% select(one_of(cov_ids)) %>% as.matrix
  ate_bal <- try(ATE::ATE(Y = this_data$y, Ti = this_data$d,
                          X = this_covm), silent = TRUE)
  if (class(ate_bal) == 'try-error') {
    data.frame(bal_omega = NA)
  } else data.frame(bal_omega = ate_bal$weights.p - ate_bal$weights.q)
}
