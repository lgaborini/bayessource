#
# Rcpp implementation wrapper
#

#' Fast computation of the marginal likelihood for the Normal-Inverted Wishart model.
#'
#' This is the numerator of the Bayes factor: assume that all observations come from the same source.
#' Implemented in C.
#'
#' See \code{\link{diwishart_inverse}} for the parametrization of the Inverted Wishart.
#' See \code{\link[bayessource]{marginalLikelihood_internal}} for further documentation.
#'
#' @param output.mcmc if TRUE output the entire chain as a {coda} object, else just return the log-ml value
#' @inheritParams marginalLikelihood_internal
#' @return the log-marginal likelihood value, or a list:
#' - `value`: the log-ml value
#' - `mcmc`: a `coda` object with the posterior samples
#' @example man-roxygen/example_marginal_likelihood.R
#' @family core functions
#' @export
#' @template gaussmv_model
#' @template InverseWishart_Press
#' @references \insertAllCited{}
#' @md
marginalLikelihood <- function(X, n.iter, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = FALSE, verbose = FALSE, Gibbs_only = FALSE) {

   # Wrap the C function
   result <- marginalLikelihood_internal(
      X = X,
      n_iter = n.iter,
      B_inv = B.inv,
      W_inv = W.inv,
      U = U,
      nw = nw,
      mu = mu,
      burn_in = burn.in,
      chain_output = output.mcmc,
      verbose = verbose,
      Gibbs_only = Gibbs_only
   )

   # Return the full MCMC output?
   if (!output.mcmc) {
      return(result$value)
   }

      # Build the coda object using the chain outputs
      # Skip the burn-in samples

      p <- ncol(X)
      # theta columns are named theta.1, ..., theta.p
      theta.mtx <- result$theta_gibbs
      colnames(theta.mtx) <- paste0("theta.", 1:p)

      # samples from W^(-1)
      W.inv.mtx <- t(apply(result$W_inv_gibbs, 3, as.numeric))
      colnames(W.inv.mtx) <- paste0("W.inv.", 1:(p^2))

      mcmc.data <- cbind(
         theta.mtx,
         W.inv.mtx
      )

      coda_object <- coda::mcmc(
         data = mcmc.data[(burn.in + 1):n.iter, ],
         start = (burn.in + 1),
         end = nrow(mcmc.data)
      )

      return(
         list(
            value = result$value,
            mcmc = coda_object
         )
      )
}



#' Fast Bayesian same source hypothesis for the Normal - Inverted Wishart model.
#'
#' Implemented in C.
#' See \code{\link[bayessource]{diwishart_inverse}} for the parametrization of the Inverted Wishart.
#' See \code{\link[bayessource]{marginalLikelihood_internal}} for further documentation.
#'
#' @param ref the reference dataset (nr * p matrix)
#' @param quest the questioned dataset (nq * p matrix)
#' @param W.inv.1 prior inverse of within-source covariance matrix (questioned items)
#' @param W.inv.2 prior inverse of within-source covariance matrix (reference items)
#' @param marginals if TRUE, also return the marginal likelihoods in the LR formula (default: FALSE)
#' @inheritParams marginalLikelihood
#' @return the log-BF value (base e), or a list with the log-BF and the computed marginal likelihoods:
#'
#' - `value`: the log-BF value (base e)
#' - `log_ml_Hp`: log-BF numerator (from reference = questioned source)
#' - `log_ml_Hd_ref`: log-BF denominator from reference source
#' - `log_ml_Hd_quest`: log-BF denominator from questioned (!= reference) source
#'
#' @export
#' @family core functions
#' @seealso marginalLikelihood
#' @template gaussmv_model
#' @template InverseWishart_Press
#' @references \insertAllCited{}
samesource_C <- function(quest, ref, n.iter, B.inv, W.inv.1, W.inv.2, U, nw, mu, burn.in, verbose = FALSE, marginals = FALSE) {

   # Wrap the C functions
   log_ml_Hp <- marginalLikelihood_internal(
      X = rbind(quest, ref),
      B_inv = B.inv,
      W_inv = W.inv.1,
      U = U, nw = nw,
      mu = mu,
      burn_in = burn.in,
      n_iter = n.iter,
      chain_output = FALSE,
      verbose = verbose
   )

   log_ml_Hd_quest <- marginalLikelihood_internal(
      X = quest,
      B_inv = B.inv,
      W_inv = W.inv.1,
      U = U, nw = nw,
      mu = mu,
      burn_in = burn.in,
      n_iter = n.iter,
      chain_output = FALSE,
      verbose = verbose
   )

   log_ml_Hd_ref <- marginalLikelihood_internal(
      X = ref,
      B_inv = B.inv,
      W_inv = W.inv.2,
      U = U, nw = nw,
      mu = mu,
      burn_in = burn.in,
      n_iter = n.iter,
      chain_output = FALSE,
      verbose = verbose
   )

   log_BF <- log_ml_Hp$value - log_ml_Hd_quest$value - log_ml_Hd_ref$value

   if (marginals) {
      return(list(
         value = log_BF,
         log_ml_Hp = log_ml_Hp$value,
         log_ml_Hd_ref = log_ml_Hd_ref$value,
         log_ml_Hd_quest = log_ml_Hd_quest$value
      ))
   } else {
      return(log_BF)
   }
}
