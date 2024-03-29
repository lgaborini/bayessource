#
# Post-processing utilities

#' Post-process Gibbs chain outputs.
#'
#' Extract and reshape MCMC samples from the posteriors for theta and W
#' .
#' Can be used to post-process outputs from \code{\link{marginalLikelihood}} when \code{output.mcmc} is TRUE.
#' Notice that it performs matrix inversions on every sample of W.
#'
#' It is also possible to compute ML estimators for the posterior means, as well as the cumulated estimates.
#' The results will be appended to the returned list.
#'
#' @param mcmc.output The coda object
#' @param compute.ML if TRUE, also compute the posterior means
#' @param cumulative if TRUE, also compute cumulative posterior means to assess precision
#' @return a named list with multidimensional arrays for W, theta, and ML estimates
#' @family core functions
#' @export
mcmc_postproc <- function(mcmc.output, compute.ML = FALSE, cumulative = TRUE) {
   stopifnot(coda::is.mcmc(mcmc.output))

   cummean <- function(x){ cumsum(x)/seq(length(x)) }

   n.samples <- nrow(mcmc.output)
   # Recover p from matrix output:
   # theta columns are named 'theta.1', 'theta.2', ..., 'theta.p'
   # count them
   p <- length(grep('^theta\\.[0-9]+$', colnames(mcmc.output)))

   theta.samples <- mcmc.output[, paste0('theta.', seq(1:p))]
   W.inv.samples <- mcmc.output[, paste0('W.inv.', seq(1:(p^2)))]
   W.inv.samples.cube <- array(W.inv.samples, dim = c(n.samples, p, p))
   W.samples.cube <- array(t(apply(W.inv.samples.cube, 1, solve)), dim = c(n.samples, p, p))

   list.out <- list(theta.samples = theta.samples, W.samples = W.samples.cube)

   if (compute.ML) {
      list.out$W.samples.ML <- apply(W.samples.cube, c(2,3), mean)
      list.out$theta.samples.ML <- as.vector(colMeans(theta.samples))

      if (cumulative) {
         list.out$W.samples.ML.cum <- apply(W.samples.cube, c(2,3), cummean)
         list.out$theta.samples.ML.cum <- apply(as.matrix(theta.samples), 2, cummean)
      }
   }
   list.out
}
