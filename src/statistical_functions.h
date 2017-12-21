#ifndef _STATISTICAL_FUNCTIONS_H
#define _STATISTICAL_FUNCTIONS_H

#include "config.h"

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

#include <limits>


// -------- Numerical functions -----------

//' Computes \eqn{log( sum_i( exp(v[i] )) )} in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double logSumExp(const arma::vec &v);

//' Computes \eqn{log( sum_i( exp(v[i] )) ) - log(n)} in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double logSumExpMean(const arma::vec &v);

//' Computes log( cumsum_i( exp(v[i]) ) ) in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::vec logCumsumExp(const arma::vec &v);

//' Computes log( cummean_i( exp(v[i]) ) ) in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::vec logCummeanExp(const arma::vec &v);

//' Upper triangular matrix inversion
//'
//' Quickly computes the inverse of a upper triangular matrix (e.g. a Cholesky factor).
//'
//' Equivalent R code:
//'
//' \code{X.chol.inv <- backsolve(r = X.chol, x = diag(p))}
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::mat inv_triangular(const arma::mat &U);

//' Compute the inverse from the upper Cholesky factor
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::mat chol2inv(const arma::mat &U_chol);

//' Upper Cholesky factor of inverse from upper Cholesky factor
//'
//' If \eqn{A = U' U}, compute \eqn{V} where \eqn{A^{(-1)} = V' V}
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::mat inv_Cholesky_from_Cholesky(const arma::mat &U);


//' log-determinant from Cholesky factor
//'
//' If \eqn{A = U' U}, compute log(det(A)) from U
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double ldet_from_Cholesky(const arma::mat &T_chol);

// -------- Statistical functions -----------

//' Generate from multivariate normal.
//'
//' Faster than \code{rmvnorm} in package \pkg{mvtnorm}. Implemented in C.
//'
//' @param n amount of samples to generate from
//' @param mu column vector for the mean
//' @param Cov covariance matrix
//' @param is_chol if TRUE, Cov is the upper Cholesky factor of Cov
//' @return a nxp matrix of samples
//' @export
// [[Rcpp::export]]
arma::mat rmvnorm(const unsigned int n,
                       const arma::colvec &mu,
                       const arma::mat &Cov,
                       const bool is_chol = false);

//' Multivariate normal density. Assumes symmetry.
//'
//' Faster than \code{dmvnorm} in package \pkg{mvtnorm}. Implemented in C.
//'
//' @param x the observation (nxp)
//' @param mean mean vector (row vector, 1xp)
//' @param Cov covariance matrix (pxp)
//' @param logd if TRUE, return the log-density
//' @param is_chol if TRUE, Cov is the upper Cholesky factor of Cov
//' @return the density in x (nx1)
//' @export
// [[Rcpp::export(rng = false)]]
arma::vec dmvnorm(const arma::mat &x,
                       const arma::rowvec &mean,
                       const arma::mat &Cov,
                       const bool logd = false,
                       const bool is_chol = false);

//' Generate random sample from Wishart (faster).
//'
//' Same code as \code{\link[stats]{rWishart}} function in package \pkg{stats}.
//'
//' @param v dof
//' @param S the scale matrix (pxp)
//' @param is_chol if TRUE, S is the upper Cholesky factor of S
//' @param return_chol if TRUE, the upper Cholesky factor is returned
//' @return a single random variate from W(v, S)
//' @export
//'
//' @template Wishart_eqn
// [[Rcpp::export]]
arma::mat rwish(const double v,
   const arma::mat &S,
   const bool is_chol = false,
   const bool return_chol = false);


//' Inverted Wishart density from the inverse (faster).
//'
//' Computes the density of an Inverted Wishart (df, Sigma) in X, by supplying (X^(-1), df, Sigma) rather than (X, df, Sigma).
//' Avoids a matrix inversion.
//'
//' Computes the pdf p_X(x) by knowing x^(-1)
//'
//' @param X_inv inverse of X (the observation)
//' @param df degrees of freedom
//' @param Sigma scale matrix
//' @param logd if TRUE, return the log-density
//' @param is_chol if TRUE, Sigma and X.inv are the upper Cholesky factors of Sigma and X.inv
//' @export
//' @template InverseWishart_Press
//' @seealso \code{\link{diwishart}}, \code{\link{dwishart}}
// [[Rcpp::export(rng = false)]]
double diwishart_inverse(const arma::mat &X_inv,
                       const double &df,
                       const arma::mat &Sigma,
                       const bool logd = false,
                       const bool is_chol = false);





#endif
