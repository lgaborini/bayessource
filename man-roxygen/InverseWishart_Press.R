#' @section Inverted Wishart parametrization (Press):
#'
#' Uses \insertCite{Press2012Applied}{bayessource} parametrization.
#'
#' \deqn{X \sim IW(\nu, S)}{X ~ IW(v, S)}
#' with \eqn{S} is a \eqn{p \times p}{p x p} matrix, \eqn{\nu > 2p}{v > 2p} (the degrees of freedom).
#'
#' Then:
#' \deqn{E[X] = \frac{S}{\nu - 2(p + 1)}}{E[X] = S/(n - 2(p + 1))}
#'
