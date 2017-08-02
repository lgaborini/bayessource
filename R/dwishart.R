
#' Wishart density.
#'
#' Density of X ~ Wishart(df, Sigma).
#'
#' @param X the observation
#' @param df degrees of freedom
#' @param Sigma scale matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, Sigma and X are the upper Cholesky factors of Sigma and X
#' @return the density in X
#' @export
#'
#' @template Wishart_eqn
dwishart <- function(X, df, Sigma, log = FALSE, is.chol = FALSE) {
   p <- ncol(X)

   # The constant
   lc0 <- -df * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df + 1 - (1:p))/2)))

   if (!is.chol) {
      X.chol <- chol(X)
      Sigma.chol <- chol(Sigma)
   } else {
      X.chol <- X
      Sigma.chol <- Sigma
   }
   ldetSigma <- 2 * sum(log(diag(Sigma.chol)))
   ldetX <- 2 * sum(log(diag(X.chol)))
   Sigma.chol.inv <- backsolve(r = Sigma.chol, x = diag(p))      # inverse of Cholesky factor (NOT chol2inv(Sigma.chol) = solve(X)!)

   # The log(trace):
   # ltrace <- sum(diag(chol2inv(Sigma.chol) %*% crossprod(X.chol)))
   # ltrace <- sum((Sigma.chol %*% solve(X.chol))^2)
   ltrace <- sum((X.chol %*% Sigma.chol.inv)^2)

   ldens <- lc0 + (df - p - 1)/2 * ldetX - (df / 2) * ldetSigma - ltrace / 2

   if (log == FALSE) { return(exp(ldens)) }
   return(ldens)
}