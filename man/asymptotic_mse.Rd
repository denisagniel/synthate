% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asymptotic_mse.R
\name{asymptotic_mse}
\alias{asymptotic_mse}
\title{Estimate asymptotic MSE}
\usage{
asymptotic_mse(V_0, V, C, deltahat)
}
\arguments{
\item{V_0}{scalar variance of theta_0}

\item{V}{covariance matrix of other candidate estimators}

\item{C}{vector of covariances of theta_0 with candidates}

\item{deltahat}{estimate of biases of candidates}
}
\value{
estimated MSE of synthetic estimator
}
\description{
Estimate asymptotic MSE
}
\examples{

k <- 5
n <- 100
theta <- c(0, runif(k-1))
Sigma <- 10*diag(k)/n
d <- theta[-1]
V_0 <- Sigma[1,1]
V <- Sigma[-1,-1]
C <- Sigma[-1,1]
asymptotic_mse(V_0, V, C, d)

}
