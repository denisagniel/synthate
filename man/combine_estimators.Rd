% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_estimators.R
\name{combine_estimators}
\alias{combine_estimators}
\title{Create a synthetic estimator by combining multiple candidate estimators.}
\usage{
combine_estimators(
  ests,
  name_0 = NULL,
  boot_ests = NULL,
  cov = NULL,
  print = FALSE,
  exclude_t0 = FALSE,
  bias_type = "raw_diff",
  ate_0 = NULL,
  n = NULL,
  ...
)
}
\arguments{
\item{ests}{one-row, p-column data frame of estimators.}

\item{name_0}{character value of the name of the presumed unbiased estimator, \eqn{\theta_0}. Default is NULL, which returns results for using each candidate estimator as \eqn{\theta_0}, one synthetic estimator for each. If \code{ate_0} is given, then it is used as \eqn{\theta_0} in place of this.}

\item{boot_ests}{p-column matrix of bootstrap estimators corresponding to the estimators in \code{ests}. If NULL, then \code{cov} must be supplied.}

\item{cov}{p x p covariance matrix of \code{ests}. Default is NULL, but if supplied \code{boot_ests} are not needed.}

\item{print}{logical indicating whether details should be printed. Default is FALSE.}

\item{exclude_t0}{logical indicating whether \eqn{\theta_0} should be considered an external estimator (not a candidate for combining with others). Default is FALSE.}

\item{bias_type}{method to compute the bias in the mean squard error. Default is \code{raw_diff}, which computes the bias as the raw difference between each of the candidate estimator and \eqn{\theta_0}. Other options to compute the bias include: \code{bootstrap} which computes the bias as the difference between the mean of the bootstrap samples and the observed value of \eqn{\theta_0}; \code{bootstrap_all} which computes the bias as the mean of the difference between the bootstrapped version of the candidate estimator and the bootstrapped version of \eqn{\theta_0}; \code{none} which assumes no bias; \code{shrunk} which computes the bias as the raw difference divided by \code{n}.}

\item{ate_0}{external value of \eqn{\theta_0}. Default is NULL, in which case \eqn{\theta_0} is taken to be \code{name_0}.}

\item{n}{sample size. Default is NULL. Needed only if \code{bias_type} is \code{shrunk}.}
}
\value{
list of three objects, including \code{ate_res} which gives results for the synthetic estimator, \code{b_res} which gives results for how the estimators were combined, and \code{C} which gives the covariance matrix of the estimators.
}
\description{
Creates a synthetic estimator by minimizing the (estimated) mean squared error of a linear combination of multiple candidate estimators.
}
\examples{

gen_mod <- generate_data(n = 100, 
                         dgp = 'ks', 
                         correct_outcome = FALSE,
                         correct_ps = TRUE)
this_data <- gen_mod$data
this_data <- estimate_scores(this_data, outcome_fm = outcome_fm,
                             ps_fm = ps_fm,
                             ps_fam = ps_fam,
                             outcome_fam = outcome_fam)
thetahat <- with(gen_mod,
                 estimate_ates(this_data,
                          ate_list,
                          cov_ids = cov_ids,
                          outcome_fm = stringr::str_c('d + ', outcome_fm),
                          outcome_fam = outcome_fam))
                          predict_delta <- function(d) {
                          gen_mod$true_ate
                          }
 resample_thetas <- with(gen_mod, resample_fn(dat = this_data,
                                dpredfn = predict_delta,
                                B = B,
                                ate_list = ate_list,
                                outcome_fm = outcome_fm,
                                ps_fm = ps_fm,
                                ps_fam = ps_fam,
                                outcome_fam = outcome_fam))
 boot_theta <- resample_thetas[[1]] \%>\% as.matrix
 
 synthetic_estimator <- combine_estimators(thetahat,
 boot_ests = boot_theta,
 name_0 = 'ate_dr',
 bias_type = 'raw_diff')
 synthetic_estimator$ate_res
}
