#' @section Inverted Wishart parametrization (Anderson):
#'
#' Uses \insertCite{Anderson2003introduction}{bayessource} parametrization.
#'
#' \deqn{X \sim IW(v, S)}{X ~ IW(v, S)}
#' with \eqn{S} is a \eqn{p \times p}{p x p} matrix, \eqn{\nu > 0}{v > 0} (the degrees of freedom).
#'
#' Then:
#'
#' \deqn{E[X] = \frac{S}{\nu - p - 1}}{E[X] = S / (v - p - 1)}
#'

