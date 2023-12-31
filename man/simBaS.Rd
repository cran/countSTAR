% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{simBaS}
\alias{simBaS}
\title{Compute Simultaneous Band Scores (SimBaS)}
\usage{
simBaS(sampFuns)
}
\arguments{
\item{sampFuns}{\code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve}
}
\value{
\code{m x 1} vector of simBaS
}
\description{
Compute simultaneous band scores (SimBaS) from Meyer et al. (2015, Biometrics).
SimBaS uses MC(MC) simulations of a function of interest to compute the minimum
alpha such that the joint credible bands at the alpha level do not include zero.
This quantity is computed for each grid point (or observation point) in the domain
of the function.
}
\note{
The input needs not be curves: the simBaS may be computed
for vectors to achieve a multiplicity adjustment.

The minimum of the returned value, \code{PsimBaS_t},
over the domain \code{t} is the Global Bayesian P-Value (GBPV) for testing
whether the function is zero everywhere.
}
