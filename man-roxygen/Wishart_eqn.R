#' @section Wishart parametrization:
#'
#' Uses \insertCite{Press2012Applied,Anderson2003introduction}{bayessource} parametrization.
#'
#' \deqn{X \sim W(\nu, S)}{X ~ W(\nu, S)}
#' with \eqn{S} is a \eqn{p \times p}{pxp} matrix, \eqn{\nu >= p}{v >= p} (the degrees of freedom).
#'
#' Then:
#' \deqn{E[X] = \nu S}{E[X] = v * S}
#'
