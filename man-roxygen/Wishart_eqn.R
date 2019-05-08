#' @section Wishart parametrization:
#'
#' Uses \insertCite{Press2012Applied,Anderson2003introduction}{bayessource} parametrization.
#'
#' \deqn{X ~ W(v, S)}
#' with \eqn{S = pxp} matrix, \eqn{v >= p} (the degrees of freedom).
#'
#' Then:
#' \deqn{E[X] = v * S}
#'
