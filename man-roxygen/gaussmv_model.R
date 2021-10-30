#' @section Normal-Inverse-Wishart model:
#'
#' Described in \insertCite{Bozza2008Probabilistic}{bayessource}.
#'
#' Observation level:
#'
#' - \deqn{X_{ij} \sim N_p(\theta_i, W_i)}{X_{ij} ~ N_p(theta_i, W_i)} (i = source, j = items from source)
#'
#' Group level:
#'
#' - \deqn{\theta_i \sim N_p(\mu, B)}{theta_i ~ N_p(\mu, B)}
#' - \deqn{W_i \sim IW_p(\nu_w, U)}{W_i ~ IW_p(nw, U)}
#'
#' Hyperparameters:
#'
#' - \deqn{B, U, \nu_w, \mu}{B, U, nw, \mu}
#'
#' Posterior samples of \eqn{\theta}{theta}, \eqn{W^{(-1)}} can be generated with a Gibbs sampler.
#'
