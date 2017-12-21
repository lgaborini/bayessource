#' @section Inverted Wishart parametrization (Anderson):
#' 
#' Uses Anderson[1] parametrization.
#' 
#' \deqn{X ~ IW(v, S)} with \eqn{S = pxp} matrix, \eqn{v > 0} (the degrees of freedom).
#' 
#' Then:
#' \deqn{E[X] = S / (v - p - 1)}
#' 
#' @references [1]T. W. Anderson, An introduction to multivariate statistical analysis, 3rd ed. Hoboken, N.J: Wiley-Interscience, 2003.
