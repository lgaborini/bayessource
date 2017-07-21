#include <limits>

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
//// [[Rcpp::plugins(cpp11)]]

// Use Cholesky factorization when possible.
// Matrices are propagated through their Cholesky factors.
#define USE_CHOLESKY false
// #define USE_CHOLESKY true

using namespace std;
using namespace Rcpp;

const double log2pi = log(2.0 * arma::datum::pi);
const double pi = arma::datum::pi;


// -------- Numerical functions -----------

//' Computes \eqn{log( sum_i( exp(v[i] )) )} in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double logSumExp(const arma::vec &v){
   double v_max = v.max();

   return (v_max + log(sum(exp(v - v_max))));
}

//' Computes \eqn{log( sum_i( exp(v[i] )) ) - log(n)} in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double logSumExpMean(const arma::vec &v){
   return (logSumExp(v) - log(v.n_elem));
}

//' Computes log( cumsum_i( exp(v[i]) ) ) in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::vec logCumsumExp(const arma::vec &v){
   double v_max = v.max();

   return (v_max + log(cumsum(exp(v - v_max))));
}

//' Computes log( cummean_i( exp(v[i]) ) ) in a stable way.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::vec logCummeanExp(const arma::vec &v){
   return (logCumsumExp(v) - log(v.n_elem));
}

//' Upper triangular matrix inversion
//'
//' Quickly computes the inverse of a upper triangular matrix (e.g. a Cholesky factor).
//'
//' Equivalent R code:
//'
//' \code{X.chol.inv <- backsolve(r = X.chol, x = diag(p))}
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::mat inv_triangular(const arma::mat &U){
   return(arma::solve(arma::trimatu(U), arma::eye(arma::size(U))));
}

//' Compute the inverse from the upper Cholesky factor
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::mat chol2inv(const arma::mat &U_chol){
   arma::mat X = solve(arma::trimatl(U_chol.t()), arma::eye(arma::size(U_chol)));
   return solve(arma::trimatu(U_chol), X);
}

//' Upper Cholesky factor of inverse from upper Cholesky factor
//'
//' If \eqn{A = U' U}, compute \eqn{V} where \eqn{A^{(-1)} = V' V} 
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
arma::mat inv_Cholesky_from_Cholesky(const arma::mat &U){
   return (arma::chol(chol2inv(U)));
}

//' log-determinant from Cholesky factor
//'
//' If \eqn{A = U' U}, compute log(det(A)) from U
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
double ldet_from_Cholesky(const arma::mat &T_chol){
   return 2 * (arma::sum(log(arma::diagvec(T_chol))));
}

// -------- Statistical functions -----------

//' Generate from multivariate normal.
//'
//' Faster than \pkg{mvtnorm}::\code{rmvnorm} (implemented in C).
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
                       const bool is_chol = false) {
   unsigned int ncols = Cov.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   arma::mat Cov_chol;

   if (is_chol)
      Cov_chol = Cov;
   else
      Cov_chol = arma::chol(Cov);

   return arma::repmat(mu, 1, n).t() + Y * Cov_chol;
}


//' Multivariate normal density. Much faster, assumes symmetry.
//'
//' Faster than \code{\link{dmvnorm_sym}}. Implemented in C.
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
                       const bool is_chol = false) {
   int n = x.n_rows;
   int p = x.n_cols;
   arma::vec out(n);

   arma::mat Cov_chol;
   if (is_chol) {
      Cov_chol = Cov;
   } else {
      Cov_chol = arma::chol(Cov);
   }

   arma::mat rooti = arma::trans(inv_triangular(Cov_chol));
   double rootisum = arma::sum(log(rooti.diag()));
   double constants = -(static_cast<double>(p)/2.0) * log2pi;

   for (int i=0; i < n; i++) {
      arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
      out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
   }

   if (logd == false) {
      out = exp(out);
   }
   return(out);
}

//' Generate random sample from Wishart (faster).
//'
//' Same code as \code{\link{rWishart}} function in package \pkg{stats}.
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
   const bool return_chol = false){

   // Dimension of returned wishart
   unsigned int p = S.n_rows;

   // The Cholesky factor
   // CC: upper triangular
   // By default, CC is upper triangular, such that CC.t()*CC = X
   arma::mat CC;
   if (is_chol) {
      CC = S;
   } else {
      CC = arma::chol(S);
   }

   // Z composition:
   // sqrt chisqs on diagonal
   // random normals above diagonal
   // 0 below diagonal
   arma::mat Z(p, p, arma::fill::zeros);

   // Fill the diagonal
   for(arma::uword i = 0; i < p; i++){
      Z.at(i,i) = sqrt(R::rchisq(v - i));
   }

   // Fill the upper matrix with random guesses
   for(arma::uword i = 0; i < p; ++i){
      for(arma::uword j = i + 1; j < p; ++j){
         Z.at(i,j) = R::rnorm(0,1);
      }
   }

   // Z is triu
   // Upper triangular * chol decomp
   // C is triu
   arma::mat C = Z * CC;

   // Return random Wishart
   if (return_chol) {
      return(C);
   } else {
      return(C.t() * C);
   }
}



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
                       const bool is_chol = false) {

   unsigned int p = Sigma.n_cols;
   double out;

   // Constants
   arma::vec vp = arma::linspace<arma::vec>(1, p, p);
   double lc0 = -(df - p - 1.) * p/2. * log(2.) - p*(p - 1.)/4. * log(pi) - arma::sum(arma::lgamma((df - p - vp)/2.));

   arma::mat X_inv_chol;
   arma::mat Sigma_chol;
   if (is_chol){
      X_inv_chol = X_inv;
      Sigma_chol = Sigma;
   } else {
      X_inv_chol = arma::chol(X_inv);
      Sigma_chol = arma::chol(Sigma);
   }

   // double lDetX = -log(arma::det(X_inv));
   // double lDetSigma = log(arma::det(Sigma));
   double lDetX = -ldet_from_Cholesky(X_inv_chol);
   double lDetSigma = ldet_from_Cholesky(Sigma_chol);

   // The trace: all equivalent expressions
   // // ltrace <- sum(diag(X.inv %*% Sigma))
   // // ltrace <- sum(X.inv * Sigma)
   // // ltrace <- sum((Sigma.chol %*% t(X.inv.chol))^2)          # how to avoid transpose?

   double ltrace = arma::trace(X_inv * Sigma);
   // double ltrace = arma::accu(X_inv % Sigma);
   // double ltrace = arma::as_scalar(arma::accu(arma::square(Sigma_chol * X_inv_chol.t())));

   out = lc0 - df * lDetX / 2. + (df - p - 1.)/2. * lDetSigma - ltrace / 2.;

   if (logd == false) {
      out = exp(out);
   }
   return(out);
}



/*** R
# Run the unit test for this script
# testthat::test_file('test_samesource_cpp.R')

*/
