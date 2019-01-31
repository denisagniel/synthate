#--------------------------------
## all data-generating mechanisms
#--------------------------------
inv.logit <- plogis

generate_data_pa <- function(n, pi = 0.2) {
  d <- rbinom(n, 1, pi)
  x1 <- rnorm(n,d*0.2)
  x2 <- rnorm(n,d*0.3)
  x3 <- rnorm(n,d*0.4)
  x4 <- rnorm(n,d*0.5)
  x5 <- rnorm(n,d*0.6)
  v1 <- rbinom(n, 1, 0.1*(1-d) + 0.168*d)
  v2 <- rbinom(n, 1, 0.2*(1-d) + 0.331*d)
  v3 <- rbinom(n, 1, 0.3*(1-d) + 0.492*d)
  v4 <- rbinom(n, 1, 0.4*(1-d) + 0.642*d)
  v5 <- rbinom(n, 1, 0.5*(1-d) + 0.776*d)
  y <- -2.8 + 1.5*x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 + 5*v1 + 4*v2 + 3*v3 + 2*v4 + 1.5*v5 + d + rnorm(n, sd = 2)
  data.frame(id = 1:n, x1, x2, x3, x4, x5, v1, v2, v3, v4, v5, y = y, d = d)
}
generate_model_pa <- function(correct_outcome, correct_ps) {
  cov_ids <- c(stringr::str_c('x', 1:5), stringr::str_c('v', 1:5))
  if (correct_outcome) {
    outcome_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  } else {
    outcome_fm <- stringr::str_c(cov_ids[-c(3)], collapse = ' + ')
  }
  if (correct_ps) {
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }  else {
    cov_ids <- c(str_c('x', c(1,2,4,5)), str_c('v', c(1,2,4,5)))
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }

  data_fn <- function(n) generate_data_pa(n, pi = 0.2)
  outcome_fam <- gaussian
  ps_fam <- binomial
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam,
       cov_ids = cov_ids, true_ate = 1)
}

generate_data_r <- function(n, pi) {
  d <- rbinom(n, 1, pi)
  x1 <- rnorm(n,d*0.2)
  x2 <- rnorm(n,d*0.3)
  x3 <- rnorm(n,d*0.4)
  x4 <- rnorm(n,d*0.5)
  x5 <- rnorm(n,d*0.6)
  v1 <- rbinom(n, 1, 0.1*(1-d) + 0.168*d)
  v2 <- rbinom(n, 1, 0.2*(1-d) + 0.331*d)
  v3 <- rbinom(n, 1, 0.3*(1-d) + 0.492*d)
  v4 <- rbinom(n, 1, 0.4*(1-d) + 0.642*d)
  v5 <- rbinom(n, 1, 0.5*(1-d) + 0.776*d)
  y <- -2.8 + 1.5*(x1 > 1) + 2*(x2 < -1) + 3*(x3 > 1)*(x1 < 1) + 4*(x4 < 1) + 5*(x5 > 0) +
    5*v1 + 4*v2 + 3*v3 + 2*v4 + 1.5*v5 + d + rnorm(n, sd = 2)
  data.frame(id = 1:n, x1, x2, x3, x4, x5, v1, v2, v3, v4, v5, y = y, d = d)
}
generate_model_r <- function(correct_outcome, correct_ps) {
  cov_ids <- c(stringr::str_c('x', 1:5), stringr::str_c('v', 1:5))
  if (correct_outcome) {
    outcome_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  } else {
    outcome_fm <- stringr::str_c(cov_ids[-c(3)], collapse = ' + ')
  }
  if (correct_ps) {
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }  else {
    cov_ids <- c(str_c('x', c(1,2,4,5)), str_c('v', c(1,2,4,5)))
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }

  data_fn <- function(n) generate_data_r(n, pi = 0.2)
  outcome_fam <- gaussian
  ps_fam <- binomial
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam)
}


generate_data_fi <- function(n, correct_ps, correct_outcome) {
  z <- rnorm(4*n) %>% matrix(n, 4)
  x <- z
  x[,1] <- exp(z[,1]/3)
  x[,2] = z[,2]/(1+exp(z[,1])) + 10
  x[,3] = z[,1]*z[,3]/25 + 0.6
  x[,4] = z[,2]+z[,4]+20
  if (correct_ps) {
    pi_0 <- cbind(1, z) %*% c(0, -1, 0.5, -0.25, -0.1)
  } else {
    pi_0 <- cbind(1, x) %*% c(0, -1, 0.5, -0.25, -0.1)
  }

  d <- rbinom(n, 1, inv.logit(pi_0))
  if (correct_outcome) {
    mu_0 <- 200 + 27.4*z[,1]*d + 13.5*z[,2] + 13.5*z[,3] + 13.5*z[,4]
  } else {
    mu_0 <- 200 + 27.4*z[,1]^2*d + 13.5*z[,2]^2 + 13.5*z[,3]^2 + 13.5*z[,4]^2
  }
  y <- rnorm(n, mu_0, 1)
  data_0 <- data.frame(id = 1:n, x = x, z = z, d = d, y = y)
  data_0
}

generate_model_fi <- function(correct_outcome, correct_ps) {
  cov_ids <- stringr::str_c('z.', 1:4)
  outcome_fm <- ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')

  data_fn <- function(n) generate_data_fi(n, correct_ps = correct_ps,
                                          correct_outcome = correct_outcome)
  outcome_fam <- gaussian
  ps_fam <- binomial
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam,
       cov_ids = cov_ids,
       true_ate = case_when(
         correct_outcome ~ 0,
         TRUE ~ 27.4
       ))
}
#
generate_data_ik <- function(n) {
  library(Matching)
  data(lalonde)
  ik_data <- dplyr::sample_n(lalonde, n, replace = TRUE)
  ik_data <- ik_data %>% mutate(
    muhat = 1 + 1.428e-4 * age^2 - 2.918e-3 * educ^2 -
      0.2275 * black - 0.8276 * hisp + 0.2071 * married -
      0.8232 * nodegr - 1.236e-9 * re74^2 +
      5.865e-10 * re75^2 - 0.04328 * u74 - 0.3804 * u75,
    pi_0 = 1 + 0.5*muhat + 0.01 * age^2 - 0.3 * educ^2 -
      0.01 * log((re74 + 0.01)^2) + 0.01 * log((re75 + 0.01)^2),
    d = rbinom(n, 1, inv.logit(pi_0)),
    l75 = log((re75 + 0.01)^2),
    ll74 = log((re74 + 0.01 + 0.7 * l75)^2),
    mu_0 = 1000*d + 0.1 * exp(0.7 * ll74),
    y = mu_0 + rnorm(n, 0, 10),
    age2 = age^2,
    educ2 = educ^2,
    l74 = log((re74 + 0.01)^2),
    re74_2 = re74^2,
    re75_2 = re75^2,
    id = 1:n
  )
}
generate_model_ik <- function(correct_outcome, correct_ps) {
  if (correct_outcome) {
    outcome_fm <- c('bs(re74)*bs(re75)')
  } else {
    outcome_fm <- 'age + educ + re74 + re75 + black +
          hisp + married + nodegr + u74 + u75'
  }
  if (correct_ps) {
    cov_ids <- c('age2', 'educ2', 'black', 'hisp', 'married',
                 'nodegr', 'u74', 'u75',
                 'l74', 'l75', 're74_2', 're75_2')
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }  else {
    cov_ids <- c('age', 'educ', 'black', 'hisp', 'married',
                 'nodegr', 'u74', 'u75',
                 're74', 're75')
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }

  data_fn <- function(n) generate_data_ik(n)
  outcome_fam <- gaussian
  ps_fam <- binomial
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam,
       cov_ids = cov_ids, true_ate = 1000)
}

generate_data_ks <- function(n, a = c(0,-1,0.5,-0.25,-0.1),
                             b = c(210,40,27.4,13.7,13.7,13.7),
                             sigma = 1) {
  z <- rnorm(4*n) %>% matrix(n, 4)
  x <- z
  x[,1] <- exp(z[,1]/2)
  x[,2] = z[,2]/(1+exp(z[,1])) + 10
  x[,3] = (z[,1]*z[,3]/25 + 0.6)^3
  x[,4] = (z[,2]+z[,4]+20)^2

  pi_0 <- cbind(1, z) %*% a
  d <- rbinom(n, 1, inv.logit(pi_0))
  mu_0 <- cbind(1, d, z) %*% b
  y <- rnorm(n, mu_0, sigma)
  data_0 <- data.frame(id = 1:n, x = x, z = z, d = d, y = y)
  data_0
}
generate_model_ks <- function(correct_outcome, correct_ps) {
  if (correct_outcome) {
    cov_ids <- stringr::str_c('z.', 1:4)
    outcome_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  } else {
    cov_ids <- stringr::str_c('x.', 1:4)
    outcome_fm <- stringr::str_c(cov_ids[-c(3)], collapse = ' + ')
  }
  if (correct_ps) {
    cov_ids <- stringr::str_c('z.', 1:4)
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }  else {
    cov_ids <- stringr::str_c('x.', 1:4)
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }

  data_fn <- function(n) generate_data_ks(n)
  outcome_fam <- gaussian
  ps_fam <- binomial
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam,
       cov_ids = cov_ids, true_ate = 40)
}

generate_data_iw <- function(n) {
  x1 <- runif(n, 1, 5)
  x3 <- x1^2
  x2 <- runif(n, 1, 5)
  x4 <- x2^2
  pi <- (1 + exp(-0.36 + 1.25*x1 + 1.25*x2 - 0.36*x1^2 - 0.35*x2^2))^-1
  t <- rbinom(n, prob= pi, size= 1)
  mu <- 2*t + exp(1 + 0.2*x1 + 0.2*x2)
  y <- rpois(n, lambda = mu)
  data.frame(id = 1:n, x.1 = x1, x.2 = x2, x.3 = x3, x.4 = x4, d = t, y = y)
}
generate_model_iw <- function(correct_outcome, correct_ps) {
  if (correct_outcome) {
    outcome_fam <- poisson
  } else {
    outcome_fam <- poisson(link = 'identity')
  }
  if (correct_ps) {
    ps_fam <- binomial
  }  else {
    ps_fam <- binomial(link = 'probit')
  }

  data_fn <- function(n) generate_data_iw(n)
  cov_ids <- stringr::str_c('x.', 1:2)
  outcome_fm <- ps_fm <- stringr::str_c(stringr::str_c('x.', 1:4), collapse = ' + ')
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam,
       cov_ids = cov_ids, true_ate = 2)
}

generate_data_ld <- function(n) {
  beta <- c(0,0.6,-0.6,0.6)
  zeta <- c(-1,1,1)
  nu <- c(0, -1, 1, -1, 2)

  tau1 <- c(1,1,-1,-1)
  tau0 <- c(-1,-1,1,1)
  Sigma <- matrix(c(1,0.5,-0.5,-0.5,
                    0.5,1,-0.5,-0.5,
                    -0.5,-0.5,1,0.5,
                    -0.5,-0.5,0.5,1), nrow=4, byrow=TRUE)
  X3 <- rbinom(n,1,0.2)
  V3 <- rbinom(n,1,0.75*X3+0.25*(1-X3))
  other_covariates <- matrix(NA, n, 4)
  other_covariates[X3 == 1,] <- mvnfast::rmvn(sum(X3 == 1), mu = tau1, sigma = Sigma)
  other_covariates[X3 == 0,] <- mvnfast::rmvn(sum(X3 == 0), mu = tau0, sigma = Sigma)
  X1 <- other_covariates[,1]
  V1 <- other_covariates[,2]
  X2 <- other_covariates[,3]
  V2 <- other_covariates[,4]

  e = inv.logit(beta[1]+beta[2]*X1+beta[3]*X2+beta[4]*X3)
  Z = rbinom(n,1,e)
  Y = nu[1] + nu[2]*X1 + nu[3]*X2 + nu[4]*X3 + nu[5]*Z +
    zeta[1]*V1 + zeta[2]*V2 + zeta[3]*V3 + rnorm(n)
  # plot(e,Y,col=ifelse(Z==1,"black","red"))
  # list(Y=Y,Z=Z,X1=X1,X2=X2,X3=X3,V1=V1,V2=V2,V3=V3)
  data.frame(id = 1:n, y = Y, d = Z, x.1 = X1, x.2 = X2, x.3 = X3,
             v.1 = V1, v.2 = V2, v.3 = V3)
}

generate_model_ld <- function(correct_outcome, correct_ps) {
  cov_ids <- c(stringr::str_c('x.', 1:3), stringr::str_c('v.', 1:3))
  if (correct_outcome) {
    outcome_fm <- ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  } else {
    outcome_fm <- stringr::str_c(cov_ids[-c(1,4)], collapse = ' + ')
  }
  if (correct_ps) {
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }  else {
    cov_ids <- c(stringr::str_c('x.', 2:3), stringr::str_c('v.', 2:3))
    ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  }

  data_fn <- function(n) generate_data_ld(n)
  outcome_fam <- gaussian
  ps_fam <- binomial
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam,
       cov_ids = cov_ids, true_ate = 2)
}

generate_data_ls <- function(n, ps_correct, outcome_correct) {
  sigma <- diag(10)
  sigma[1,5] <- sigma[5,1] <- sigma[3,8] <- sigma[8,3] <- 0.2
  sigma[2,6] <- sigma[6,2] <- sigma[3,9] <- sigma[9,3] <- 0.9
  w <- mvnfast::rmvn(n, mu = rep(0, 10), sigma = sigma)
  w[,c(1,3,5,6,8,9)] <- ifelse(w[,c(1,3,5,6,8,9)] > 0, 1, 0)
  w_ps <- cbind(1, w, w[,2]^2, w[,4]^2, w[,7]^2, w[,1]*w[,3],
                w[,2]*w[,4], w[,3]*w[,5], w[,6]*w[,4],
                w[,5]*w[,7], w[,1]*w[,6], w[,2]*w[,3],
                w[,3]*w[,4], w[,5]*w[,4], w[,5]*w[,6])
  alpha <- c(-1.897, 0.8, -0.25, 0.6, -0.4, -0.8, -0.5, 0.7, 0, 0, 0,
               -0.25, -0.4, 0.7, 0.5*0.8,
               0.7*(-0.25), 0.5*0.6, 0.7*(-0.4),
               0.5*(-0.8), 0.5*0.8, 0.7*(-0.25),
               0.5*0.6, 0.5*(-0.4), 0.5*(-0.8))
  w_o <- cbind(1, w, w[,2]^2, w[,4]^2, w[,10]^2,
               w[,1]*w[,3],
               w[,2]*w[,4], w[,3]*w[,8], w[,9]*w[,4],
               w[,8]*w[,10], w[,1]*w[,9], w[,2]*w[,3],
               w[,3]*w[,4], w[,8]*w[,4], w[,8]*w[,9])
  beta <- c(-1.386, 0.3, -0.36, -0.73, -0.2, 0, 0, 0, 0.71, -0.19, 0.26,
            -0.36, -0.2, 0.26,
            0.5*0.3, 0.7*(-0.36), 0.5*(-0.73), 0.7*(-0.2),
            0.5*0.71, 0.5*0.3, 0.7*(-0.36),
            0.5*(-0.73), 0.5*(-0.2), 0.5*(0.71))
  pi_0 <- w_ps %*% alpha
  d <- rbinom(n, 1, plogis(pi_0))
  y <- w_o %*% beta - 0.4*d + rnorm(n, sd = 0.1)
  data.frame(id = 1:n, y = y, d = d, ww = w, w_ps = w_ps, w_o = w_o)
}

generate_model_ls <- function(correct_outcome, correct_ps) {

  if (correct_outcome) {
    outcome_fm <- stringr::str_c(stringr::str_c('w_o.', 2:24),
                                 collapse = ' + ')
  } else {
    outcome_fm <- stringr::str_c(stringr::str_c('ww.', 1:10),
                                 collapse = ' + ')
  }
  if (correct_ps) {
    cov_ids <- stringr::str_c('w_ps.', 2:24)
  }  else {
    cov_ids <- stringr::str_c('ww.', 1:10)
  }
  ps_fm <- stringr::str_c(cov_ids, collapse = ' + ')
  data_fn <- function(n) generate_data_ls(n)
  outcome_fam <- gaussian
  ps_fam <- binomial
  list(data_fn = data_fn, outcome_fm = outcome_fm,
       outcome_fam = outcome_fam,
       ps_fm = ps_fm, ps_fam = ps_fam,
       cov_ids = cov_ids, true_ate = -0.4)
}

generate_data_sj <- function(n, gamma_c = 0.5, compliance_p = 0.5, compliance_effect = 0.5, alpha_c = 1, alpha_n = 0.5, lambda_c = 1, lambda_n = 1, gamma_n = 0, sigma_c = 1, sigma_n = 1) {
  # browser()
  x <- rnorm(n) ## covariate
  z <- rbinom(n, prob = 0.5, size = 1) ## randomization indicator
  if (compliance_effect == 0) {
    pi_c <- plogis(qlogis(compliance_p)) ## probability of compliance
  } else {
    pi_c <- plogis(log(compliance_effect)*x + qlogis(compliance_p)) ## probability of compliance
  }

  c <- rbinom(n, prob = pi_c, size = 1) ## compliance indicator
  epsilon <- rnorm(n)


  y <- alpha_n + (alpha_c - alpha_n)*c + gamma_n*z + (gamma_c - gamma_n)*c*z +
    lambda_n*x + (lambda_c - lambda_n)*c*x + sigma_n*epsilon + (sigma_c - sigma_n)*c*epsilon

  sim_data <- data.frame(id = 1:n,
                         x = x,
                         z = z,
                         c = c,
                         s = 1*(z == 1 & c == 1), ## receipt of treatment
                         y = y)
  sim_data
}

generate_data <- function(n, dgp, correct_outcome = TRUE, correct_ps = TRUE) {
  genfn <- get(paste0('generate_model_', dgp))
  mod <- genfn(correct_outcome, correct_ps)
  data <- mod$data_fn(n)
  list(data = data,
       dgp = dgp,
       outcome_fm = mod$outcome_fm,
       outcome_fam = mod$outcome_fam,
       ps_fm = mod$ps_fm,
       ps_fam = mod$ps_fam,
       cov_ids = mod$cov_ids,
       true_ate = mod$true_ate)
}
