#include "statistical_functions.h"


using namespace std;
using namespace Rcpp;

const double log2pi = log(2.0 * arma::datum::pi);
const double pi = arma::datum::pi;


// -------- Numerical functions -----------


double logSumExp(const arma::vec &v){
   double v_max = v.max();

   return (v_max + log(sum(exp(v - v_max))));
}


double logSumExpMean(const arma::vec &v){
   return (logSumExp(v) - log(v.n_elem));
}


arma::vec logCumsumExp(const arma::vec &v){
   double v_max = v.max();

   return (v_max + log(cumsum(exp(v - v_max))));
}


arma::vec logCummeanExp(const arma::vec &v){
   return (logCumsumExp(v) - log(v.n_elem));
}


arma::mat inv_triangular(const arma::mat &U){
   return(arma::solve(arma::trimatu(U), arma::eye(arma::size(U))));
}

arma::mat inv_sympd_tol(const arma::mat &U_sympd){
   return arma::inv_sympd(arma::symmatu(U_sympd));
}

arma::mat chol2inv(const arma::mat &U_chol){
   arma::mat X = solve(arma::trimatl(U_chol.t()), arma::eye(arma::size(U_chol)));
   return solve(arma::trimatu(U_chol), X);
}


arma::mat inv_Cholesky_from_Cholesky(const arma::mat &U){
   return (arma::chol(chol2inv(U)));
}


double ldet_from_Cholesky(const arma::mat &T_chol){
   return 2 * (arma::sum(log(arma::diagvec(T_chol))));
}

// -------- Statistical functions -----------


arma::mat rmvnorm(const unsigned int n,
                       const arma::colvec &mu,
                       const arma::mat &Cov,
                       const bool is_chol) {
   unsigned int ncols = Cov.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   arma::mat Cov_chol;

   if (is_chol)
      Cov_chol = Cov;
   else
      Cov_chol = arma::chol(Cov);

   return arma::repmat(mu, 1, n).t() + Y * Cov_chol;
}



arma::vec dmvnorm(const arma::mat &x,
                       const arma::rowvec &mean,
                       const arma::mat &Cov,
                       const bool logd,
                       const bool is_chol) {
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


arma::mat rwish(const double v,
   const arma::mat &S,
   const bool is_chol,
   const bool return_chol){
   // Sample from Wishart
   // Using Anderson/Press parametrisation.

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




double diwishart_inverse(const arma::mat &X_inv,
                       const double &df,
                       const arma::mat &Sigma,
                       const bool logd,
                       const bool is_chol) {

   unsigned int p = Sigma.n_cols;
   double out;

   const double df_min = 2*(p + 1) + 1;
   if (df < df_min) {
      Rcpp::stop("Error: Inverted Wishart is degenerate (%u df, dimension %u, minimum df %u).", df, p, df_min);
   }

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
