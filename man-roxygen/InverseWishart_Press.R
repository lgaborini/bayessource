#' @section Inverted Wishart parametrization (Press):
#' 
#' Uses Press[Press2012] parametrization.
#' 
#' \deqn{X ~ IW(v, S)}
#' with \eqn{S = pxp} matrix, \eqn{v > 2p} (the degrees of freedom).
#' 
#' Then:
#' \deqn{E[X] = S / (v - 2(p + 1))}
#' 
#' @references [Press2012] S. J. Press, Applied multivariate analysis: using Bayesian and frequentist methods of inference. Courier Corporation, 2012.
