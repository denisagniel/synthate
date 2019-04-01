prematch <- function(this_data, match_fns, fn_nms, M) {
  this_data <- this_data %>% mutate(index = 1:nrow(this_data))
  match_list <- list()
  for (i in 1:length(match_fns)) {
    match_list[[i]] <- prematch_(this_data, match_fns[[i]], fn_nms[[i]], M[i])
  }
  left_join(this_data, Reduce(full_join, match_list))
}

prematch_ <- function(this_data, match_fn, fn_nm, M) {
  # browser()
  this_match <- match_fn(this_data)
  K_ds <- get_K(this_match, fn_nm, M)
  K_ds
}

get_K <- function(match_fit, nm, M) {
  # browser()
  ctrl_tab <- table(match_fit$index.control)
  case_tab <- table(match_fit$index.treated)

  k_ds <- tibble(index = as.numeric(c(names(ctrl_tab),
                           names(case_tab))),
                 K = c(ctrl_tab,
                       case_tab)/M
  )
  colnames(k_ds)[2] <- paste0(nm, '_K')
  k_ds
}


generic_match <- function(this_data, matching_vars,
                              M = 5, nm, ...) {
  match <- with(this_data, Matching::Match(Y = y, Tr = d, X = matching_vars,
                                           estimand = 'ATE', version = 'fast',
                                           ties = FALSE, M = M, ...))
  match
}

match_ps <- function(this_data, ...) {
  generic_match(this_data, this_data$ps, nm = 'match_ps')
}
match_prog <- function(this_data, ...) {
  generic_match(this_data, this_data$prog_score, nm = 'match_prog')
}
match_both <- function(this_data, ...) {
  generic_match(this_data, cbind(this_data$ps, this_data$prog_score), nm = 'match_both')
}

caliper_ps <- function(this_data, ...) {
  generic_match(this_data, this_data$ps, nm = 'cal_match_ps',
                    caliper = 0.05)
}
caliper_prog <- function(this_data, ...) {
  generic_match(this_data, this_data$prog_score, nm = 'cal_match_prog',
                    caliper = 0.05)
}
caliper_both <- function(this_data, ...) {
  generic_match(this_data, cbind(this_data$ps, this_data$prog_score),
                    nm = 'cal_match_both', caliper = 0.2)
}

caliper_nn <- function(this_data, ...) {
  generic_match(this_data, this_data$ps, nm = 'cal_match_nn',
                    caliper = 0.05, replace = FALSE, M = 1)
}

caliper_logit <- function(this_data, ...) {
  generic_match(this_data, log(this_data$ps/(1-this_data$ps)), nm = 'cal_match_logit',
                    caliper = 0.05, replace = FALSE, M = 1)
}
