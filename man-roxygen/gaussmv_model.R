#' @section Multivariate Gaussian model:
#'
#' Observation level:
#' \deqn{X_{ij} ~ N(theta, W)}
#'
#' Group level:
#' \deqn{theta ~ N(\mu, B)}
#' \deqn{W ~ IW(nw, U)}
#'
#' Hyperparameters:
#' \eqn{B, U, nw, \mu}
#'
#' A Gibbs sampler supplies theta, W^{(-1)}.
