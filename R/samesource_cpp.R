#
# Rcpp implementation wrapper
#

library(coda)

#' Fast Bayesian same source hypothesis. Gaussian MV.
#' Implemented in C (faster)
#'
#' @param dati the dataset
#' @param n.iter number of MC iterations
#' @param B.inv prior inverse of between covariance matrix
#' @param W.inv prior inverse of within covariance matrix
#' @param U covariance matrix for the mean
#' @param nw degrees of freedom
#' @param mu prior mean
#' @param burn.in burn-in iterations
#' @param output.mcmc output the entire chain
#' @param verbose if TRUE, also output all posterior samples
#'
#' @return the LR value, or a list(LR value, posterior samples)
#' @export
#' @template gaussmv_model
#'
samesource_C <- function(dati, n.iter, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = FALSE, verbose = FALSE) {

   result <- samesource_C_internal(dati, n.iter, B.inv, W.inv, U, nw, mu, burn.in, chain_output = output.mcmc, verbose = verbose)

   if (output.mcmc) {

      theta.mtx <- result$theta_gibbs
      colnames(theta.mtx) <- paste0('theta.', 1:p)

      W.inv.mtx <- t(apply(result$W_inv_gibbs, 3, as.numeric))
      colnames(W.inv.mtx) <- paste0('W.inv.', 1:(p^2))

      mcmc.data <- cbind(
         theta.mtx,
         W.inv.mtx)

      return(
         list(
            LR.num = result$LR.num,
            mcmc = mcmc(data = mcmc.data[(burn.in + 1):n.iter, ], start = (burn.in + 1), end = nrow(mcmc.data))
      ))
   } else {
      return(result$LR.num)
   }
}
