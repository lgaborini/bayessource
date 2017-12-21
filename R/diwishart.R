#' Inverted Wishart density, parametrization according to Press.
#' 
#' Implemented in R (slow).
#'
#' @param X observation
#' @param df degrees of freedom
#' @param Sigma scale matrix
#' @param log if TRUE, return the log-density
#' @param is.chol if TRUE, Sigma and X are the upper Cholesky factors of Sigma and X
#' @return the density in X
#' @export
#' @seealso \code{\link{diwishart_inverse}}, \code{\link{dwishart}}
#' @rdname statistical_functions
#'
#' @template InverseWishart_Press
diwishart <- function(X, df, Sigma, log = FALSE, is.chol = FALSE) {
   p <- ncol(X)

   # The constant
   lc0 <- -(df - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((df - p - (1:p))/2)))

   if (!is.chol) {
      X.chol <- chol(X)
      Sigma.chol <- chol(Sigma)
   } else {
      X.chol <- X
      Sigma.chol <- Sigma
   }
   # ldetSigma <- log(det(Sigma))
   # ldetX <- log(det(X))
   ldetSigma <- 2 * sum(log(diag(Sigma.chol)))
   ldetX <- 2 * sum(log(diag(X.chol)))

   # Inverse of Cholesky factor
   # NOT inverse from Cholesky: chol2inv(X.chol) = solve(X)!
   # NOT Cholesky of inverse
   X.chol.inv <- backsolve(r = X.chol, x = diag(p))

   # The trace: all equivalent expressions
   # ltrace <- sum(diag(solve(X) %*% Sigma))
   # ltrace <- sum(solve(X) * Sigma)
   # ltrace <- sum(diag(chol2inv(X.chol) %*% crossprod(Sigma.chol)))
   # ltrace <- sum((Sigma.chol %*% solve(X.chol))^2)
   ltrace <- sum((Sigma.chol %*% X.chol.inv)^2)

   ldens <- lc0 + (df - p - 1)/2 * ldetSigma - (df / 2) * ldetX - ltrace / 2

   if (log == FALSE) { return(exp(ldens)) }
   return(ldens)
}

