# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Check whether Cholesky factorization is used.
#' @keywords internal
#' @family C++ functions
isCholeskyOn <- function() {
    .Call('_bayessource_isCholeskyOn', PACKAGE = 'bayessource')
}

#' Fast computation of the Normal-IW marginal likelihood.
#'
#' This is the numerator of the Bayes factor: assume that all observations come from the same source.
#' To be called by the R wrapper.
#'
#' @param X the observation matrix (\eqn{n \times p}{(n x p)}: n = observation, p = variables)
#'  @param U covariance matrix for the mean (\eqn{p \times p}{p x p})
#' @param n_iter number of MCMC iterations excluding burn-in
#' @param burn_in number of MCMC burn-in iterations
#' @param B_inv prior inverse of between-source covariance matrix
#' @param W_inv initialization for prior inverse of within-source covariance matrix
#' @param nw degrees of freedom
#' @param mu prior mean (\eqn{p \times 1}{p x 1})
#' @param chain_output if true, output the entire chain as a list (ML-value, samples from theta, samples from W_inv)
#' @param verbose if TRUE, be verbose
#' @param Gibbs_only if TRUE, only return the Gibbs posterior samples. Implies `chain_output = TRUE`.
#'
#' @template gaussmv_model
#' @keywords internal
#' @family C++ functions
#' @family core functions
marginalLikelihood_internal <- function(X, n_iter, B_inv, W_inv, U, nw, mu, burn_in, chain_output = FALSE, verbose = FALSE, Gibbs_only = FALSE) {
    .Call('_bayessource_marginalLikelihood_internal', PACKAGE = 'bayessource', X, n_iter, B_inv, W_inv, U, nw, mu, burn_in, chain_output, verbose, Gibbs_only)
}

#' Computes \eqn{log( sum_i(exp(v_i)) )} in a stable way.
#'
#' @keywords internal
#' @family C++ functions
#' @family math functions
logSumExp <- function(v) {
    .Call('_bayessource_logSumExp', PACKAGE = 'bayessource', v)
}

#' Computes \eqn{log( sum_i( exp(v_i)) ) - log(n)} in a stable way.
#'
#' @keywords internal
#' @family C++ functions
#' @family math functions
logSumExpMean <- function(v) {
    .Call('_bayessource_logSumExpMean', PACKAGE = 'bayessource', v)
}

#' Computes log( cumsum_i( exp(v_i) )) in a stable way.
#'
#' @keywords internal
#' @family C++ functions
#' @family math functions
logCumsumExp <- function(v) {
    .Call('_bayessource_logCumsumExp', PACKAGE = 'bayessource', v)
}

#' Computes log( cummean_i( exp(v_i) )) in a stable way.
#'
#' @keywords internal
#' @family C++ functions
#' @family math functions
logCummeanExp <- function(v) {
    .Call('_bayessource_logCummeanExp', PACKAGE = 'bayessource', v)
}

#' Upper triangular matrix inversion
#'
#' Quickly computes the inverse of a upper triangular matrix (e.g. a Cholesky factor).
#'
#' Equivalent R code:
#'
#' \code{X.chol.inv <- backsolve(r = X.chol, x = diag(p))}
#' @keywords internal
#' @family C++ functions
#' @family math functions
inv_triangular <- function(U) {
    .Call('_bayessource_inv_triangular', PACKAGE = 'bayessource', U)
}

#' Compute the inverse of a symmetric positive definite matrix
#'
#' Compute the inverse of a symmetric positive definite matrix.
#' Does not output warnings on default symmetry tolerance: basically, symmetry is forced.
#'
#' @references https://github.com/RcppCore/RcppArmadillo/issues/257
#' @family C++ functions
#' @keywords internal
#' @family math functions
inv_sympd_tol <- function(U_sympd) {
    .Call('_bayessource_inv_sympd_tol', PACKAGE = 'bayessource', U_sympd)
}

#' Compute the inverse from the upper Cholesky factor
#'
#' @family C++ functions
#' @keywords internal
#' @family math functions
chol2inv <- function(U_chol) {
    .Call('_bayessource_chol2inv', PACKAGE = 'bayessource', U_chol)
}

#' Upper Cholesky factor of inverse from upper Cholesky factor
#'
#' If \eqn{A = U' U}, compute \eqn{V} where \eqn{A^{(-1)} = V' V}
#'
#' @family C++ functions
#' @keywords internal
#' @family math functions
inv_Cholesky_from_Cholesky <- function(U) {
    .Call('_bayessource_inv_Cholesky_from_Cholesky', PACKAGE = 'bayessource', U)
}

#' log-determinant from Cholesky factor
#'
#' If \eqn{A = U' U}, compute \eqn{\log{\det{A}}}{log(det(A))} from U
#'
#' @keywords internal
#' @family C++ functions
#' @family math functions
ldet_from_Cholesky <- function(T_chol) {
    .Call('_bayessource_ldet_from_Cholesky', PACKAGE = 'bayessource', T_chol)
}

#' Generate from multivariate normal.
#'
#' Faster than [mvtnorm::rmvnorm()]. Implemented in C.
#'
#' @param n amount of samples to generate from
#' @param mu column vector for the mean
#' @param Cov covariance matrix
#' @param is_chol if TRUE, Cov is the upper Cholesky factor of Cov
#' @return a nxp matrix of samples
#' @family C++ functions
#' @family statistical functions
#' @export
rmvnorm <- function(n, mu, Cov, is_chol = FALSE) {
    .Call('_bayessource_rmvnorm', PACKAGE = 'bayessource', n, mu, Cov, is_chol)
}

#' Multivariate normal density. Assumes symmetry.
#'
#' Faster than [mvtnorm::dmvnorm()]. Implemented in C.
#'
#' @param x the observation (nxp matrix)
#' @param mean mean vector (row vector, 1xp)
#' @param Cov covariance matrix (pxp)
#' @param logd if TRUE, return the log-density
#' @param is_chol if TRUE, Cov is the upper Cholesky factor of Cov
#' @return the density in x (nx1)
#' @family C++ functions
#' @family statistical functions
#' @export
dmvnorm <- function(x, mean, Cov, logd = FALSE, is_chol = FALSE) {
    .Call('_bayessource_dmvnorm', PACKAGE = 'bayessource', x, mean, Cov, logd, is_chol)
}

#' Generate random sample from Wishart (faster).
#'
#' Same code as \code{\link[stats]{rWishart}} function in package \pkg{stats}.
#'
#' @param v dof
#' @param S the scale matrix (pxp)
#' @param is_chol if TRUE, S is the upper Cholesky factor of S
#' @param return_chol if TRUE, the upper Cholesky factor is returned
#' @return a single random variate from W(v, S)
#' @family C++ functions
#' @family statistical functions
#' @family Wishart functions
#' @export
#' @template Wishart_eqn
#' @references \insertAllCited{}
rwish <- function(v, S, is_chol = FALSE, return_chol = FALSE) {
    .Call('_bayessource_rwish', PACKAGE = 'bayessource', v, S, is_chol, return_chol)
}

#' Inverted Wishart density from the inverse (faster).
#'
#' Computes the pdf p_X(x) by knowing x^(-1)
#'
#' Computes the density of an Inverted Wishart (df, Sigma) in x, by supplying (x^(-1), df, Sigma) rather than (x, df, Sigma).
#' Avoids a matrix inversion.
#'
#' @param X_inv inverse of X (the data)
#' @param df degrees of freedom of the Inverted Wishart
#' @param Sigma scale matrix of the Inverted Wishart
#' @param logd if TRUE, return the log-density
#' @param is_chol if TRUE, Sigma and X_inv are the upper Cholesky factors of Sigma and X_inv
#' @family C++ functions
#' @family statistical functions
#' @family Wishart functions
#' @export
#' @template InverseWishart_Press
#' @seealso \code{\link{diwishart}}, \code{\link{dwishart}}
#' @references \insertAllCited{}
diwishart_inverse <- function(X_inv, df, Sigma, logd = FALSE, is_chol = FALSE) {
    .Call('_bayessource_diwishart_inverse', PACKAGE = 'bayessource', X_inv, df, Sigma, logd, is_chol)
}

