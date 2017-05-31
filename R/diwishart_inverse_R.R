
#' Inverted Wishart density from the inverse
#'
#' Computes the density of an Inverted Wishart (df, Sigma) in X, by supplying (X^(-1), df, Sigma) rather than (X, df, Sigma).
#' Avoids a matrix inversion.
#'
#' Computes the pdf p_X(x) by knowing x^(-1)
#'
#' @param X.inv inverse of X (the observation)
#' @param df degrees of freedom
#' @param Sigma scale matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, Sigma and X.inv are the upper Cholesky factors of Sigma and X.inv
#' @return the density in X
#' @keywords internal
#'
#' @template InverseWishart_Press
diwishart_inverse_R <- function(X.inv, df, Sigma, log = FALSE, is.chol = FALSE) {
   p <- ncol(X.inv)

   # The constant
   lc0 <- -(df - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df - p - (1:p))/2)))

   if (!is.chol) {
      X.inv.chol <- chol(X.inv)
      Sigma.chol <- chol(Sigma)
   } else {
      X.inv.chol <- X.inv
      Sigma.chol <- Sigma
   }

   # Cholesky factor of X: inverse-transpose of Cholesky factor of X.inv
   # X.chol <- t(backsolve(r = X.inv.chol, x = diag(p)))

   # ldetSigma <- log(det(Sigma))
   # ldetX <- log(det(X))
   ldetSigma <- 2 * sum(log(diag(Sigma.chol)))
   ldetX.inv <- 2 * sum(log(diag(X.inv.chol)))
   ldetX <- -ldetX.inv

   # The trace:
   # ltrace <- sum(diag(X.inv %*% Sigma))
   # ltrace <- sum(X.inv * Sigma)
   ltrace <- sum((Sigma.chol %*% t(X.inv.chol))^2)          # how to avoid transpose?

   ldens <- lc0 + (df - p - 1)/2 * ldetSigma - (df / 2) * ldetX - ltrace / 2

   if (log == FALSE) { return(exp(ldens)) }
   return(ldens)
}
