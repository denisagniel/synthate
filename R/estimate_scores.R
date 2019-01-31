estimate_scores <- function(this_data, outcome_fm, outcome_fam = gaussian,
                            ps_fm, ps_fam = binomial, wt = NULL, ...) {
  # browser()
  ps_fm <- as.formula(stringr::str_c('d ~ ', ps_fm))
  prog_fm <- as.formula(stringr::str_c('y ~ ', outcome_fm))
  n <- nrow(this_data)

  if (!is.null(wt)) {
    browser()
    this_data <- this_data %>% mutate(w = wt)
    ps_design <- survey::svydesign(ids = ~1, weights = ~w, data = this_data, formula = ps_fm)
    ps_fit <- survey::svyglm(ps_fm, data = this_data, family = ps_fam, design = ps_design)
    prog_design <- survey::svydesign(ids = ~1,
                                     weights = ~w,
                                     data = this_data %>% filter(d == 0),
                                     formula = prog_fm)
    lm_ini <- lm(prog_fm, data = this_data %>% filter(d == 0),
                 weights = w)
    if (!identical(outcome_fam, gaussian)) {
      prog_fit <- survey::svyglm(prog_fm, data = this_data %>% filter(d == 0),
                      family = outcome_fam, start = coef(lm_ini), design = prog_design)
    } else prog_fit <- lm_ini
  } else {
    this_data <- this_data %>% mutate(w = 1)
    ps_fit <- glm(ps_fm,
                data = this_data, family=ps_fam)
    lm_ini <- lm(prog_fm, data = this_data %>% filter(d == 0))
    if (!identical(outcome_fam, gaussian)) {
      prog_fit <- glm(prog_fm, data = this_data %>% filter(d == 0),
                      family = outcome_fam, start = coef(lm_ini), etastart = predict(lm_ini))
    } else prog_fit <- lm_ini
  }

  this_data %>%
    mutate(ps = predict(ps_fit, type = 'response',
                        newdata = this_data),
           pi = ifelse(d == 1, ps, 1 - ps),
           ps_strat = as.numeric(cut(ps, 5)),
           prog_score = predict(prog_fit, type = 'response',
                                newdata = this_data))
}
