#include "config.h"
#include "statistical_functions.h"
#include "samesource.h"




using namespace std;
using namespace Rcpp;



//' Fast Bayesian same source hypothesis. Gaussian MV.
//' To be called by the R wrapper.
//'
//' @param X the dataset
//' @param n_iter number of MC iterations
//' @param B_inv prior inverse of between covariance matrix
//' @param W_inv prior inverse of within covariance matrix
//' @param U covariance matrix for the mean
//' @param nw degrees of freedom
//' @param mu prior mean
//' @param burn_in burn-in iterations
//' @param chain_output output the entire chain
//' @param verbose if TRUE, be verbose
//' @param Gibbs_only if TRUE, only return the Gibbs posterior samples. Implies chain_output = TRUE.
//' @param n_cores if > 1, use OpenMP on the specified amount of cores (experimental).
//'
//' @template gaussmv_model
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List marginalLikelihood_internal(
      const arma::mat &X,
      const unsigned int n_iter,
      const arma::mat &B_inv,
      const arma::mat &W_inv,
      const arma::mat &U,
      const double nw,
      const arma::vec &mu,
      const unsigned int burn_in,
      const bool chain_output = false,
      const bool verbose = false,
      const bool Gibbs_only = false,
      const unsigned int n_cores = 1){


   // try/catch block
   BEGIN_RCPP

   // Dimensions
   unsigned int nr = X.n_rows;
   unsigned int p = X.n_cols;

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

   // OpenMP
   omp_set_num_threads(n_cores);
   if (verbose){
      Rcout << "Using OpenMP with " << n_cores << " threads." << endl;
   }

   // Compute non-normalized sample covariance
   arma::rowvec bary = arma::mean(X, 0);
   arma::mat S(p,p);
   S.zeros();
   for (arma::uword j = 0; j < nr; ++j){
      S += (X.row(j) - bary).t() * (X.row(j) - bary);
   }

   // Initialize Gibbs chain
   // chain over theta_g (mu_upd_g, B_upd_g), W_inv_g (U_upd_g)
   //
   // We initialize W_inv_g from the hyperprior parameter W_inv
   // B, mu are updated using the posterior formulae
   // theta is sampled using its definition
   // As a consequence, only W_inv_g is fully arbitrary
   W_inv_g = W_inv;

   if (verbose){
      Rcout << "Gibbs chain started." << endl;
      Rcout << "Parameters: " << endl;
      Rcout << "  using Cholesky? " << USE_CHOLESKY << endl;
   }
   int print_mod = 100;

   for (arma::uword i = 0; i < n_iter; ++i){
      // if (verbose && (i % print_mod == 0)) Rcout << ".";

      B_upd_g = arma::inv_sympd(B_inv + nr * W_inv_g);
      mu_upd_g = B_upd_g * (B_inv * mu + nr * W_inv_g * bary.t());

      if (USE_CHOLESKY){
         B_upd_g_chol = arma::chol(B_upd_g);
         theta_g = rmvnorm(1, mu_upd_g, B_upd_g_chol, true);
      } else {
         theta_g = rmvnorm(1, mu_upd_g, B_upd_g, false);
      }


      // Cholesky: has no effect here
      U_upd_g = nr * (theta_g - bary).t() * (theta_g - bary) + U + S;

      // Sample W from Inverted Wishart (nwstar, U_upd_g) = sample from Wishart (nwstar, solve(U_upd_g)), then invert
      if (USE_CHOLESKY){
         U_upd_g_chol = arma::chol(U_upd_g);
         W_inv_g_chol = rwish(nwstar, inv_Cholesky_from_Cholesky(U_upd_g_chol), true, true);
      } else {
         // Using stats::rWishart (verified, same result!)
         // W_inv_g = Rcpp::as<arma::cube>(rWishart(1, nwstar, arma::inv_sympd(U_upd_g))).slice(0);
         W_inv_g = rwish(nwstar, arma::inv_sympd(U_upd_g), false, false);
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
         logf = arma::sum(dmvnorm(X, theta_g, inv_Cholesky_from_Cholesky(W_inv_g_chol), true, true));
      } else {
         logf = arma::sum(dmvnorm(X, theta_g, arma::inv_sympd(W_inv_g), true, false));
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
      Rcout << "Maximum likelihood converged to (\\theta^*, \\W^*): logf* = " << logf_star << endl;
   }

   if (Gibbs_only){
      // Convert Cholesky to full
      if (USE_CHOLESKY){
         for (arma::uword i = 0; i < n_iter; ++i) {
            W_inv_gibbs.slice(i) = W_inv_gibbs.slice(i).t() * W_inv_gibbs.slice(i);
         }
      }

      return(List::create(
            _["LR.num"] = NumericVector::create(NA_REAL),
            _["W_inv_gibbs"] = W_inv_gibbs,
            _["theta_gibbs"] = theta_gibbs
      ));
   }

   // Estimate the posterior ordinates ================

   if (verbose) Rcout << "Computing posterior ordinates." << endl;

   double lpihat_theta_star;
   double lpihat_W_star;

   #pragma omp parallel for schedule(static) private(lpihat_theta_star, lpihat_W_star)
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
         lpihat_theta_star = arma::as_scalar(dmvnorm(theta_star, arma::conv_to< arma::rowvec >::from(mu_upd_g), B_upd_g_chol, true, true));
         lpihat_W_star = diwishart_inverse(W_inv_star_chol, nwstar, U_upd_g_chol, true, true);
      } else {
         lpihat_theta_star = arma::as_scalar(dmvnorm(theta_star, arma::conv_to< arma::rowvec >::from(mu_upd_g), B_upd_g, true, false));
         lpihat_W_star = diwishart_inverse(W_inv_star, nwstar, U_upd_g, true, false);
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
   lpihat_theta_star = logSumExpMean(lpihat_theta_star_samples.subvec(burn_in + 1, n_iter - 1));
   lpihat_W_star = logSumExpMean(lpihat_W_star_samples.subvec(burn_in + 1, n_iter - 1));

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
      lpi_theta_star = arma::as_scalar(dmvnorm(theta_star, arma::conv_to< arma::rowvec >::from(mu), arma::chol(B), true, true));
   } else {
      lpi_theta_star = arma::as_scalar(dmvnorm(theta_star, arma::conv_to< arma::rowvec >::from(mu), B, true, false));
   }

   double lpi_W_star;
   if (USE_CHOLESKY) {
      lpi_W_star = diwishart_inverse(W_inv_star_chol, nw, arma::chol(U), true, true);
   } else {
      lpi_W_star = diwishart_inverse(W_inv_star, nw, U, true, false);
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
   if (USE_CHOLESKY && chain_output) {
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


   return(List::create(_["LR.num"] = 1));

   // try/catch block
   END_RCPP
}



/*** R
# Run the unit test for this script
# testthat::test_file('test_samesource_cpp.R')
*/
