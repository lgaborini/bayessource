#
# Rcpp implementation wrapper
#

#' Fast Bayesian marginal likelihood for the Normal - Inverted Wishart model.
#'
#' Implemented in C.
#' See \code{diwishart_inverse} for the parametrization of the Inverted Wishart.
#'
#' @param X the dataset
#' @param n.iter number of MC iterations
#' @param B.inv prior inverse of between covariance matrix
#' @param W.inv prior inverse of within covariance matrix
#' @param U covariance matrix for the mean
#' @param nw degrees of freedom
#' @param mu prior mean
#' @param burn.in burn-in iterations
#' @param output.mcmc output the entire chain
#' @param verbose if TRUE, be verbose
#' @param Gibbs_only other parameters passed to \code{marginalLikelihood_internal}
#'
#' @return the log-marginal likelihood value, or a list(the log-ml value, coda object with the posterior samples)
#' @export
#' @template gaussmv_model
#' @template InverseWishart_Press
#'
marginalLikelihood <- function(X, n.iter, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = FALSE, verbose = FALSE, Gibbs_only = FALSE) {

   # Wrap the C function
   result <- marginalLikelihood_internal(X, n.iter, B.inv, W.inv, U, nw, mu, burn.in, chain_output = output.mcmc, verbose = verbose, Gibbs_only = Gibbs_only)

   if (output.mcmc) {
      # Build the coda object using the chain outputs
      # Skip the burn-in samples

      p <- ncol(X)
      # theta columns are named theta.1, ..., theta.p
      theta.mtx <- result$theta_gibbs
      colnames(theta.mtx) <- paste0('theta.', 1:p)

      # samples from W^(-1)
      W.inv.mtx <- t(apply(result$W_inv_gibbs, 3, as.numeric))
      colnames(W.inv.mtx) <- paste0('W.inv.', 1:(p^2))

      mcmc.data <- cbind(
         theta.mtx,
         W.inv.mtx)

      return(
         list(
            LR.num = result$LR.num,
            mcmc = coda::mcmc(data = mcmc.data[(burn.in + 1):n.iter, ], start = (burn.in + 1), end = nrow(mcmc.data))
      ))
   } else {
      return(result$LR.num)
   }
}



#' Fast Bayesian same source hypothesis for the Normal - Inverted Wishart model.
#'
#' Implemented in C.
#' See \code{diwishart_inverse} for the parametrization of the Inverted Wishart.
#'
#' @param quest the questioned dataset
#' @param ref the reference dataset
#' @param n.iter number of MC iterations
#' @param B.inv prior inverse of between covariance matrix
#' @param W.inv.1 prior inverse of within covariance matrix (questioned items)
#' @param W.inv.2 prior inverse of within covariance matrix (reference items)
#' @param U covariance matrix for the mean
#' @param nw degrees of freedom
#' @param mu prior mean
#' @param burn.in burn-in iterations
#' @param verbose if TRUE, be verbose
#'
#' @return the log-LR value
#' @export
#' @template gaussmv_model
#' @template InverseWishart_Press
#'
samesource_C <- function(quest, ref, n.iter, B.inv, W.inv.1, W.inv.2, U, nw, mu, burn.in, verbose = FALSE) {
   # Wrap the C functions
   LR.num <- marginalLikelihood_internal(rbind(quest, ref), n.iter, B.inv, W.inv.1, U, nw, mu, burn.in, chain_output = FALSE, verbose = verbose)
   LR.den.1 <- marginalLikelihood_internal(quest, n.iter, B.inv, W.inv.1, U, nw, mu, burn.in, chain_output = FALSE, verbose = verbose)
   LR.den.2 <- marginalLikelihood_internal(ref, n.iter, B.inv, W.inv.2, U, nw, mu, burn.in, chain_output = FALSE, verbose = verbose)

   LR.num$LR.num - LR.den.1$LR.num - LR.den.2$LR.num
}
