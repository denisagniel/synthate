random_grouping <- function(this_data, K_random = 20, ate_list, outcome_fm, ps_fm, ps_fam = binomial, outcome_fam = gaussian, cov_ids) {
  # cut treatment/control separately
  index_treat = which(this_data$d==1)
  index_ctr = which(this_data$d==0)
  treat_random = cut(sample(length(index_treat),length(index_treat),replace=F),breaks=K_random,labels=FALSE)
  ctr_random = cut(sample(length(index_ctr),length(index_ctr),replace=F),breaks=K_random,labels=FALSE)

  r_theta = matrix(NA , K_random , length(thetahat))
  colnames(r_theta) <- colnames(thetahat)

  for (k_random in 1:K_random){
    random_data <- estimate_scores(this_data[c(index_treat[treat_random==k_random],
                                               index_ctr[ctr_random==k_random]),],
                                   outcome_fm = outcome_fm,
                                   ps_fm = ps_fm,
                                   ps_fam = ps_fam,
                                   outcome_fam = outcome_fam)

    r_theta[k_random,] = estimate_ates(random_data, ate_list,
                                       cov_ids = cov_ids,
                                       outcome_fm = stringr::str_c('d + ', outcome_fm),
                                       outcome_fam = outcome_fam) %>% unlist
  }
  # scale for number of random subsets
  r_theta/sqrt(K_random)
}

