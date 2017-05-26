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

// Computes log( \sum_i( exp(v[i] )) ) in a stable way.
// [[Rcpp::export(.logSumExp_arma, rng = false)]]
double logSumExp_arma(const arma::vec &v){
   double v_max = v.max();

   return (v_max + log(sum(exp(v - v_max))));
}

// Computes log( \sum_i( exp(v[i] )) ) - log(n) in a stable way.
// [[Rcpp::export(.logSumExpMean_arma, rng = false)]]
double logSumExpMean_arma(const arma::vec &v){
   return (logSumExp_arma(v) - log(v.n_elem));
}

// Computes log( \cumsum_i( exp(v[i] )) ) in a stable way.
// [[Rcpp::export(.logCumsumExp_arma, rng = false)]]
arma::vec logCumsumExp_arma(const arma::vec &v){
   double v_max = v.max();

   return (v_max + log(cumsum(exp(v - v_max))));
}

// Computes log( \cummean_i( exp(v[i] )) ) in a stable way.
// [[Rcpp::export(.logCumsumExpmean_arma, rng = false)]]
arma::vec logCumsumExpmean_arma(const arma::vec &v){
   return (logCumsumExp_arma(v) - log(v.n_elem));
}

// Upper triangular matrix inversion
//
// R: X.chol.inv <- backsolve(r = X.chol, x = diag(p))
//
// [[Rcpp::export(.inv_triangular, rng = false)]]
arma::mat inv_triangular(const arma::mat &U){
   return(arma::solve(arma::trimatu(U), arma::eye(arma::size(U))));
}

// Compute the inverse from the upper Cholesky factor
//
// [[Rcpp::export(.chol2inv, rng = false)]]
arma::mat chol2inv(const arma::mat &U_chol){
   arma::mat X = solve(arma::trimatl(U_chol.t()), arma::eye(arma::size(U_chol)));
   return solve(arma::trimatu(U_chol), X);
}

// Cholesky factor of inverse from Cholesky factor
//
// [[Rcpp::export(.inv_Cholesky_from_Cholesky, rng = false)]]
arma::mat inv_Cholesky_from_Cholesky(const arma::mat &U){
   return (arma::chol(chol2inv(U)));
}

// log-determinant from Cholesky factor
// [[Rcpp::export(rng = false)]]
double ldet_from_Cholesky(const arma::mat &T_chol){
   return 2 * (arma::sum(log(arma::diagvec(T_chol))));
}

// -------- Statistical functions -----------

// Generate from multivariate normal
// [[Rcpp::export(rmvnorm_arma)]]
arma::mat rmvnorm_arma(const unsigned int n,
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


// Density of multivariate normal
// [[Rcpp::export(dmvnorm_arma, rng = false)]]
arma::vec dmvnorm_arma(const arma::mat &x,
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


// Generate from Wishart
//
// Anderson/Press parametrization
// As R::rWishart function.
//
// Only the upper Cholesky factor of the scale matrix can be considered.
// Only the upper Cholesky factor of the sampled Wishart can be returned.
//
// [[Rcpp::export(rwish_arma)]]
arma::mat rwish_arma(const double v,
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



// Density of Inverse Wishart from inverse
// X ~ IWishart(X, df, Sigma)
// Using Press parametrization
//
// Computes the pdf p_X(x) by knowing x^(-1)
//
// [[Rcpp::export(diwishart_inverse_arma, rng = false)]]
double diwishart_inverse_arma(const arma::mat &X_inv,
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




// [[Rcpp::export(name = .samesource_C)]]
Rcpp::List samesource_C(
      const arma::mat &dati,
      const unsigned int n_iter,
      const arma::mat &B_inv,
      const arma::mat &W_inv,
      const arma::mat &U,
      const double nw,
      const arma::vec &mu,
      const unsigned int burn_in,
      const bool chain_output = false,
      const bool verbose = false) {

   // Dimensions
   unsigned int nr = dati.n_rows;
   unsigned int p = dati.n_cols;


   // GIBBS SAMPLER

   // Sample from Wishart using stats::rWishart
   // Function rWishart("rWishart", "package:stats");

   // Posterior dof
   double nwstar = nr + nw;

   // Priors
   arma::mat B = arma::inv_sympd(B_inv);

   // Gibbs updates
   arma::mat B_upd_g;
   arma::vec mu_upd_g;
   arma::mat U_upd_g;

   // Gibbs random samples
   arma::mat W_inv_g;
   arma::rowvec theta_g;

   // Cholesky factors
   arma::mat B_upd_g_chol;
   arma::mat U_upd_g_chol;
   arma::mat W_inv_g_chol;

   // Outputs Gibbs chain
   // Chain for \theta
   arma::cube B_upd_gibbs(p, p, n_iter);
   arma::mat mu_upd_gibbs(p, n_iter);
   arma::mat theta_gibbs(n_iter, p);                  // r.v.
   // Chain for W
   arma::cube U_upd_gibbs(p, p, n_iter);
   arma::cube W_inv_gibbs(p, p, n_iter);              // r.v.

   // Output likelihoods: store to compute mean of log likelihoods
   arma::vec lpihat_theta_star_samples(n_iter);
   arma::vec lpihat_W_star_samples(n_iter);

   // Maximum likelihood
   double logf;
   double logf_star = -arma::datum::inf;
   arma::rowvec theta_star;
   arma::mat W_inv_star;
   // The Cholesky factor
   arma::mat W_inv_star_chol;

   // Compute non-normalized sample covariance
   arma::rowvec bary = arma::mean(dati, 0);
   arma::mat S(p,p);
   S.zeros();
   for (arma::uword j = 0; j < nr; ++j){
      S += (dati.row(j) - bary).t() * (dati.row(j) - bary);
   }

   // Initialize Gibbs chain
   // chain over theta_g (mu_upd_g, B_upd_g), W_inv_g (U_upd_g)
   W_inv_g = W_inv;

   if (verbose){
      Rcout << "Gibbs chain started." << endl;
      Rcout << "Parameters: " << endl;
      Rcout << "  USE_CHOLESKY: " << USE_CHOLESKY << endl;
   }
   int print_mod = 100;

   for (arma::uword i = 0; i < n_iter; ++i){
      // if (verbose && (i % print_mod == 0)) Rcout << ".";

      B_upd_g = arma::inv_sympd(B_inv + nr * W_inv_g);
      mu_upd_g = B_upd_g * (B_inv * mu + nr * W_inv_g * bary.t());

      if (USE_CHOLESKY){
         B_upd_g_chol = arma::chol(B_upd_g);
         theta_g = rmvnorm_arma(1, mu_upd_g, B_upd_g_chol, true);
      } else {
         theta_g = rmvnorm_arma(1, mu_upd_g, B_upd_g, false);
      }

      U_upd_g = nr * (theta_g - bary).t() * (theta_g - bary) + U + S;

      // Sample W from Inverted Wishart (nwstar, U_upd_g) = sample from Wishart (nwstar, solve(U_upd_g)), then invert
      if (USE_CHOLESKY){
         U_upd_g_chol = arma::chol(U_upd_g);
         W_inv_g_chol = rwish_arma(nwstar, inv_Cholesky_from_Cholesky(U_upd_g_chol), true, true);
      } else {
         // Using stats::rWishart (verified, same result!)
         // W_inv_g = Rcpp::as<arma::cube>(rWishart(1, nwstar, arma::inv_sympd(U_upd_g))).slice(0);
         W_inv_g = rwish_arma(nwstar, arma::inv_sympd(U_upd_g), false, false);
      }



      // Save Gibbs outputs
      mu_upd_gibbs.col(i) = mu_upd_g;
      theta_gibbs.row(i) = theta_g;                // r.v.
      if (USE_CHOLESKY){
         // B
         B_upd_gibbs.slice(i) = B_upd_g_chol;
         // U, W
         U_upd_gibbs.slice(i) = U_upd_g_chol;
         W_inv_gibbs.slice(i) = W_inv_g_chol;              // r.v.
      } else {
         // B
         B_upd_gibbs.slice(i) = B_upd_g;
         // U, W
         U_upd_gibbs.slice(i) = U_upd_g;
         W_inv_gibbs.slice(i) = W_inv_g;              // r.v.
      }


      // Compute the likelihood
      // ... automatically
      if (USE_CHOLESKY){
         logf = arma::sum(dmvnorm_arma(dati, theta_g, inv_Cholesky_from_Cholesky(W_inv_g_chol), true, true));
      } else {
         logf = arma::sum(dmvnorm_arma(dati, theta_g, arma::inv_sympd(W_inv_g), true, false));
      }

      if (logf > logf_star){
         // if (verbose)
         //    Rcout << "[iter " << i << "] New maximum: " << logf << ", previous: " << logf_star << ", delta: " << logf - logf_star << endl;
         theta_star = theta_g;
         logf_star = logf;
         if (USE_CHOLESKY){
            W_inv_star_chol = W_inv_g_chol;
         } else {
            W_inv_star = W_inv_g;
         }
      }

      // Update priors with posteriors
      // U_g = U_upd_g;
      // B_inv_g = B_upd_inv_g;
      // mu_g = mu_upd_g;

   }

   if (verbose){
      Rcout << "Maximum likelihood converged to (\\theta^*, \\W^*): logf = " << logf << endl;
   }

   // Estimate the posterior ordinates ================

   if (verbose) Rcout << "Computing posterior ordinates." << endl;

   double lpihat_theta_star;
   double lpihat_W_star;

   for (unsigned int i = burn_in + 1; i < n_iter; ++i){

      // Load posterior hyperparameters from Gibbs chain
      mu_upd_g = mu_upd_gibbs.col(i);

      if (USE_CHOLESKY){
         B_upd_g_chol = B_upd_gibbs.slice(i);
         U_upd_g_chol = U_upd_gibbs.slice(i);
      } else {
         B_upd_g = B_upd_gibbs.slice(i);
         U_upd_g = U_upd_gibbs.slice(i);
      }

      if (USE_CHOLESKY){
         lpihat_theta_star = arma::as_scalar(dmvnorm_arma(theta_star, arma::conv_to< arma::rowvec >::from(mu_upd_g), B_upd_g_chol, true, true));
         lpihat_W_star = diwishart_inverse_arma(W_inv_star_chol, nwstar, U_upd_g_chol, true, true);
      } else {
         lpihat_theta_star = arma::as_scalar(dmvnorm_arma(theta_star, arma::conv_to< arma::rowvec >::from(mu_upd_g), B_upd_g, true, false));
         lpihat_W_star = diwishart_inverse_arma(W_inv_star, nwstar, U_upd_g, true, false);
      }

      // logsumexp: save the samples
      // notice: 0-based indexing, the last position is n_iter - 1
      lpihat_theta_star_samples[i] = lpihat_theta_star;
      lpihat_W_star_samples[i] = lpihat_W_star;
   }
   if (verbose) Rcout << "Computing posterior ordinate constants." << endl;

   // The posterior marginals as MC mean
   // lpihat_theta_star = lpihat_theta_star/(n_iter - burn_in);
   // lpihat_W_star = lpihat_W_star/(n_iter - burn_in);
   // lpihat_theta_star = arma::sum(lpihat_theta_star_samples)/(n_iter - burn_in);
   // lpihat_W_star = arma::sum(lpihat_W_star_samples)/(n_iter - burn_in);
   // The correct log-means:
   lpihat_theta_star = logSumExpMean_arma(lpihat_theta_star_samples.subvec(burn_in + 1, n_iter - 1));
   lpihat_W_star = logSumExpMean_arma(lpihat_W_star_samples.subvec(burn_in + 1, n_iter - 1));

   // Joint posterior on \theta, W: independent parameters
   double lpihat_psi_star = lpihat_theta_star + lpihat_W_star;


   if (verbose){
      Rcout << "  lpihat_theta_star: " << lpihat_theta_star << endl;
      Rcout << "  lpihat_W_star: " << lpihat_W_star << endl;
   }



// ##########################################
// ##### COMPUTE THE PRIOR ORDINATES ########
// ##########################################

   if (verbose){
      Rcout << "Computing prior ordinates." << endl;
   }

   double lpi_theta_star;
   if (USE_CHOLESKY) {
      lpi_theta_star = arma::as_scalar(dmvnorm_arma(theta_star, arma::conv_to< arma::rowvec >::from(mu), arma::chol(B), true, true));
   } else {
      lpi_theta_star = arma::as_scalar(dmvnorm_arma(theta_star, arma::conv_to< arma::rowvec >::from(mu), B, true, false));
   }

   double lpi_W_star;
   if (USE_CHOLESKY) {
      lpi_W_star = diwishart_inverse_arma(W_inv_star_chol, nw, arma::chol(U), true, true);
   } else {
      lpi_W_star = diwishart_inverse_arma(W_inv_star, nw, U, true, false);
   }
   double lpi_psi_star = lpi_theta_star + lpi_W_star;

   if (verbose){
      Rcout << "  lpi_theta_star: " << lpi_theta_star << endl;
      Rcout << "  lpi_W_star: " << lpi_W_star << endl;
   }


// ##########################################
// ##### OBTAIN THE MARGINAL DENSITY ########
// ##########################################

   double lmhatHp_num = logf_star + lpi_psi_star - lpihat_psi_star;      // Equation (8)


   if (verbose){
      Rcout << "Computing marginal density." << endl;
      Rcout << "  logf_star: " << logf_star << endl;
      Rcout << "  lpi_psi_star: " << lpi_psi_star << endl;
      Rcout << "  lpihat_psi_star: " << lpihat_psi_star << endl;
      Rcout << "LR (num):  lmhatHp_num: " << lmhatHp_num << endl;
   }


   // return(List::create(
   //       _["B_upd_gibbs"] = B_upd_gibbs,
   //       _["mu_upd_gibbs"] = mu_upd_gibbs,
   //       _["U_upd_gibbs"] = U_upd_gibbs,
   //       _["W_inv_gibbs"] = W_inv_gibbs,
   //       _["theta_gibbs"] = theta_gibbs,
   //       _["theta_star"] = theta_star,
   //       _["W_inv_star"] = W_inv_star,
   //       _["logf_star"] = logf_star,
   //       _["lpihat_theta_star"] = lpihat_theta_star,
   //       _["lDetW_inv_star"] = lDetW_inv_star,
   //       _["lc0star"] = lc0star,
   //       _["lpihat_W_star"] = lpihat_W_star,
   //       _["lpihat_psi_star"] = lpihat_psi_star,
   //       _["lpi_theta_star"] = lpi_theta_star,
   //       _["lpi_W_star"] = lpi_W_star,
   //       _["lpi_psi_star"] = lpi_psi_star,
   //       _["lmhatHp_num"] = lmhatHp_num
   // ));

   // Return the full matrices for the chain instead of Cholesky factors
   if (USE_CHOLESKY) {
      for (arma::uword i = 0; i < n_iter; ++i) {
         W_inv_gibbs.slice(i) = W_inv_gibbs.slice(i).t() * W_inv_gibbs.slice(i);
      }
   }

   if (verbose) Rcout << "Returning." << endl;

   if (chain_output) {
      return(List::create(
            _["LR.num"] = lmhatHp_num,
            _["W_inv_gibbs"] = W_inv_gibbs,
            _["theta_gibbs"] = theta_gibbs
      ));
   } else {
      return(List::create(_["LR.num"] = lmhatHp_num));
   }


}



/*** R
# Run the unit test for this script
testthat::test_file('test_samesource_cpp.R')

*/
