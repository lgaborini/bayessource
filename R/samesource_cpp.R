#
# Rcpp implementation wrapper
#

library(coda)

# Rcpp::sourceCpp('samesource.cpp', embeddedR = FALSE)         # source and skip tests
# Rcpp::sourceCpp('samesource.cpp')         # source and run tests

samesource_C <- function(dati, n.iter, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = FALSE, verbose = FALSE) {

   result <- .samesource_C(dati, n.iter, B.inv, W.inv, U, nw, mu, burn.in, chain_output = output.mcmc, verbose = verbose)

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
