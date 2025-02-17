negloglik <- function(pars, pars_baseline=NULL, obs, dist="pois", ...) {

  inc <- incidence(pars=pars, pars_baseline=pars_baseline, ...)$inc
  # multiple age groups lead to multiple columns
  # even if inc has one column, rowSums leads to a vector
  inc = rowSums(inc)

  if (length(inc) < length(obs)) {
    stop("Model predictions must be >= observations")
  }

  model = inc[(length(inc)-length(obs)+1):length(inc)]
# extract the incidence that will be evaluated against data
  # note that the early part of the outbreak are not detected in the model
  if (length(model) != length(obs)) {
    stop("The length differs between the model and observation")
  }
  if (any(model < 0)) {
    ll <- -Inf
  }
  else {
    if (dist == "pois") {
      ll <- sum(dpois(obs, model, log=TRUE), na.rm=T)
    } else { # negative binomial distribution
      # ll <- sum(dnbinom(obs, size=pars[5], mu=model, log=TRUE))
    }
  }
  return(-ll)
}
