remove_effect_from_data <- function(data, est, population = 'all') {
  if (population == 'all') {
    w1 <- 1/2
    w2 <- 1/2
  } else if (population == 'treated') {
    w1 <- 1
    w2 <- 0
  } else if (population == 'control') {
    w1 <- 0
    w2 <- 1
  } else stop('Population must be one of "all", "treated", or "control".')
    data %>%
      mutate(y_old = y,
             y = case_when(d == 1 ~ y - w1*est,
                           d == 0 ~ y + w2*est))
}
