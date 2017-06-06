#' @section Normal-Inverse-Wishart model[1]:
#'
#' Observation level:
#' \deqn{X_{ij} ~ N(theta_i, W_i)} (i = population, j = items in population)
#'
#' Group level:
#' \deqn{theta_i ~ N(\mu, B)}
#' \deqn{W_i ~ IW(nw, U)}
#'
#' Hyperparameters:
#' \eqn{B, U, nw, \mu}
#'
#' A Gibbs sampler supplies theta, W^{(-1)}.
#' 
#' @references [1]S. Bozza, F. Taroni, R. Marquis, and M. Schmittbuhl, “Probabilistic evaluation of handwriting evidence: likelihood ratio for authorship,” Journal of the Royal Statistical Society: Series C (Applied Statistics), vol. 57, no. 3, pp. 329–341, Jun. 2008.
