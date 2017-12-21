#' @section Inverted Wishart parametrization (Press):
#' 
#' Uses Press[1] parametrization.
#' 
#' \deqn{X ~ IW(v, S)}
#' with \eqn{S = pxp} matrix, \eqn{v > 2p} (the degrees of freedom).
#' 
#' Then:
#' \deqn{E[X] = S / (v - 2(p + 1))}
#' 
#' @references [1]S. J. Press, Applied multivariate analysis: using Bayesian and frequentist methods of inference. Courier Corporation, 2012.
