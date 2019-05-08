#' @section Inverted Wishart parametrization (Press):
#'
#' Uses \insertCite{Press2012Applied}{bayessource} parametrization.
#'
#' \deqn{X ~ IW(v, S)}
#' with \eqn{S = pxp} matrix, \eqn{v > 2p} (the degrees of freedom).
#'
#' Then:
#' \deqn{E[X] = S / (v - 2(p + 1))}
#'
