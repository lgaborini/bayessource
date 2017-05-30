#' @section Multivariate Gaussian model:
#'
#' Observation level:
#' \deqn{X_{ij} \sim N(theta, W)}
#'
#' Group level:
#' \deqn{theta \sim N(\mu, B)}
#' \deqn{W \sim IW(nw, U)}
#'
#' Hyperparameters:
#' \eqn{B, U, nw, \mu}
#'
#' A Gibbs sampler supplies theta, W^{(-1)}.
