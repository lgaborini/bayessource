#ifndef _STATISTICAL_FUNCTIONS_H
#define _STATISTICAL_FUNCTIONS_H


// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>



// -------- Numerical functions -----------

double logSumExp_C(const arma::vec &v);

double logSumExpMean_C(const arma::vec &v);

arma::vec logCumsumExp_C(const arma::vec &v);

arma::vec logCumsumExpmean_C(const arma::vec &v);

arma::mat inv_triangular(const arma::mat &U);

arma::mat chol2inv(const arma::mat &U_chol);

arma::mat inv_Cholesky_from_Cholesky(const arma::mat &U);

double ldet_from_Cholesky(const arma::mat &T_chol);

// -------- Statistical functions -----------

arma::mat rmvnorm_C(const unsigned int n,
                       const arma::colvec &mu,
                       const arma::mat &Cov,
                       const bool is_chol = false);


arma::vec dmvnorm_C(const arma::mat &x,
                       const arma::rowvec &mean,
                       const arma::mat &Cov,
                       const bool logd = false,
                       const bool is_chol = false);

arma::mat rwish_C(const double v,
   const arma::mat &S,
   const bool is_chol = false,
   const bool return_chol = false);


double diwishart_inverse_C(const arma::mat &X_inv,
                       const double &df,
                       const arma::mat &Sigma,
                       const bool logd = false,
                       const bool is_chol = false);





#endif
