library(rrtm)
library(BayesianTools)

obs_raw <- read.csv("my_spectra.csv")
obs <- obs_raw[obs_raw$wavelength %in% 400:2500, "reflectance"]
loglike <- function(params) {
  mod <- rrtm::prospect4(params[1], params[2], params[3], params[4])
  sum(dnorm(obs, mod, params[5], log = TRUE))
}

prior_density <- function(params) {
  dunif(params[1], 1, 10, log = TRUE) + dnorm(params[2], 40, 15, log = TRUE) + ...
}
prior_sampler <- function() {
  c(runif(1, 1, 10), rnorm(1, 40, 15), ...)
}
prior <- BayesianTools::createPrior(density = prior_density, sampler = prior_sampler, lower = c(1, 0, 0, 0, 0))
setup <- BayesianTools::createBayesianSetup(likelihood = loglike, prior = prior)
fit <- BayesianTools::runMCMC(setup)


