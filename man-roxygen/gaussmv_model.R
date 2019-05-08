#' @section Normal-Inverse-Wishart model\insertCite{Bozza2008Probabilistic}{bayessource}:
#'
#' Observation level:
#' \deqn{X_{ij} ~ N_p(theta_i, W_i)} (i = population, j = items in population)
#'
#' Group level:
#' \deqn{theta_i ~ N_p(\mu, B)}
#' \deqn{W_i ~ IW_p(nw, U)}
#'
#' Hyperparameters:
#' \eqn{B, U, nw, \mu}
#'
#' A Gibbs sampler supplies theta, W^{(-1)}.
#'
