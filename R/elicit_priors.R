#' Elicit priors and initialization from background dataset
#'
#'
#' Here we fix the hyperparameters for priors on \eqn{theta_i} and \eqn{W_i}, i.e., \eqn{B}, \eqn{U}, \eqn{mu} and \eqn{n_w}.
#'
#' An appropriate initialization value for \eqn{W^{-1}_i}, \eqn{i=1,2} is also generated.
#'
#' Notice that we have three chains in the LR computations, the same initialization is used thrice.
#'
#' ## Priors
#'
#' `use.priors`:
#' - `'ML'`: maximum likelihood estimation as in \insertCite{Bozza2008Probabilistic}{bayessource}
#' - `'vague'`: low-information priors
#'
#'   `U` is `alpha*diag(p)`, B is `beta*diag(p)`, `mu` is `mu0`.
#'
#'   By default `alpha = 1, beta = 1, mu0 = 0`
#'
#' The Wishart dofs `nw` are set as small as possible without losing full rank:
#' \eqn{nw = 2*(p + 1) + 1}
#'
#' ## Initialization
#'
#' Generate values for \eqn{W^{-1}_1}, \eqn{W^{-1}_2} that will be used as starting values for the Gibbs chains.
#' Two methods:
#'
#' `use.init`:
#'
#' - `'random'`: initialize according to the model
#'
#'    + generate one \eqn{W_i \sim IW_p(nw, U)}{W_i ~ IW_p(nw, U)}
#'    + then invert: \eqn{W^{-1}_i = (W_i)^{-1}}
#'
#' - `'vague'`: low-information initialization
#'
#'   + \eqn{W_i = alpha_init + beta_init * diag(p)}
#'
#' ## Parameters
#'
#' Some constants can be changed by passing the new values to `...` :
#'
#' - `use.priors = 'vague'`: accepts parameters
#'
#'    + `alpha` (default: `1`)
#'    + `beta` (default: `1`)
#'    + `mu0` (default: `0`)
#'
#' - `use.init = 'vague'`: accepts parameters
#'
#'    + `alpha_init` (default: `1`)
#'    + `beta_init` (default: `100`)
#'
#' ## Return values
#'
#' A list of variables:
#'
#' - `mu`: the global mean
#' - `B.inv`: the between prior covariance, as the inverse
#' - `U`: the Inverted Wishart scale matrix
#' - `nw`: the Inverted Wishart dof
#' - `W.inv.1`, `W.inv.2`: the within cov initializations, as inverses
#'
#' @export
#' @param df.background the background dataset
#' @param col.variables columns with variables
#' @param col.item column with item id
#' @param use.priors see details
#' @param use.init see details
#' @param ... additional variables for priors, init (see the description)
#' @return a list of variables
#' @family core functions
#' @example man-roxygen/example_make_priors_and_init.R
#' @md
make_priors_and_init <- function(
   df.background,
   col.variables,
   col.item,
   use.priors = c('ML', 'vague'),
   use.init = c('random', 'vague'),
   ...) {

   use.priors <- match.arg(use.priors)
   use.init <- match.arg(use.init)

   p <- length(col.variables)

   # Capture current dots: optional parameters and preset constants
   dots <- list(...)

   # minimum dof s.t. E[W] is defined, with W ~ Wishart (or Inverse Wishart)
   nw.min <- 2*(p + 1) + 1

   # Perform maximum likelihood estimation on the background dataset
   if (use.priors == 'ML') {
      # On background dataset
      # obtain B, mu, W
      WB <- two.level.multivariate.calculate.UC(df.background, col.variables, col.item)

      mu <- t(WB$all.means)
      B.inv <- chol2inv(chol(WB$B))
      nw <- nw.min

      # Choose U from W, nw
      #
      # W is generated as follows:
      # W ~ IW(nw, U)
      # We have W0 (sample estimate for W from the background)
      #
      # Set U prior s.t.
      # E[W] = U / (nw - 2*p - 2) = W0        # Press

      U <- WB$W * (nw - 2*p - 2)           # mean of Inverse Wishart (Press)

   }

   # Vague priors
   if (use.priors == 'vague') {

      dots.default <- list(alpha = 1, beta = 1, mu0 = 0)
      # Modified optional parameters
      dots.now <- utils::modifyList(dots.default, dots)

      # Generic priors
      # U is the prior on the within-source covariance matrix
      U <- dots.now$alpha * diag(p)
      mu <- rep(dots.now$mu0, p)
      # B.inv is the prior between-sources precision matrix
      B.inv <- (1 / dots.now$beta) * diag(p)

      nw <- nw.min
   }


   # End of priors
   priors <- list()
   priors$mu <- mu
   priors$U <- U
   priors$B.inv <- B.inv
   priors$nw <- nw

   # ## Chain initialization
   #
   # Here we initialize the chain by supplying $W_1^{-1}$ and $W_2^{-1}$.
   # The $\theta_i$ are generated accordingly to the model.
   #
   # Three methods:
   #
   # - sample from the hierarchical model using the hyperparameters $n_w$, $U$
   # - initialize from fixed values
   # - cheat and use the known values

   # Initialize according to the model
   # Warning: RNG is used here!
   if (use.init == 'random') {

      # W.1 <- riwish(nw.exact, U)
      W.1 <- riwish_Press(nw, U)
      W.inv.1 <- chol2inv(chol(W.1))

      # is there a better option?
      W.inv.2 <- W.inv.1
      rm(W.1)
   }

   # Initialize at a fixed value
   if (use.init == 'vague') {

      dots.default <- list(alpha_init = 1, beta_init = 100)
      dots.now <- utils::modifyList(dots.default, dots)

      # W.1 is a covariance matrix, not a precision!
      W.1 <- dots.now$alpha_init*matrix(1, p, p) + dots.now$beta_init*diag(p)
      W.inv.1 <- chol2inv(chol(W.1))
      W.inv.2 <- W.inv.1
      rm(W.1)
   }

   # Save the initialization
   priors$W.inv.1 <- W.inv.1
   priors$W.inv.2 <- W.inv.2

   priors
}

#' Get minimum degrees of freedom for Inverted Wishart
#'
#' Returns minimum value for degrees of freedom such that the Inverted Wishart has a well-defined expected value.
#'
#' Uses \insertCite{Press2012Applied}{bayessource} parametrization.
#'
#' \deqn{X \sim IW(\nu, S)}{X ~ IW(nw, S)}
#'
#' with \eqn{S = pxp} matrix, \eqn{nw > 2p} (the degrees of freedom).
#'
#' Then:
#'
#' \deqn{E[X] = S / (\nu - 2(p + 1))}{E[X] = S / (nw - 2(p + 1))}
#'
#' Finally, the minimum dof \eqn{\nu}{nw} are: \eqn{\nu = 2(p + 1) + 1}{nw = 2*(p + 1) + 1}.
#'
#' @param p number of variables
#' @return minimum dof
#' @family core functions
#' @family Wishart functions
#' @export
#' @references \insertAllCited{}
get_minimum_nw_IW <- function(p) {
   stopifnot(is.numeric(p))
   if (p < 1) { stop('p must be >= 1')}

   nw.min <- 2*(p + 1) + 1
   nw.min
}



#' Generate random sample from Inverted Wishart.
#'
#' Using \insertCite{Press2012Applied}{bayessource} parametrization.
#'
#' @param v dof (\eqn{> 2p})
#' @param S the scale matrix (pxp)
#' @return a single random variate from IW(v, S)
#' @template InverseWishart_Press
#' @family R functions
#' @family statistical functions
#' @family Wishart functions
#' @export
#' @references \insertAllCited{}
riwish_Press <- function(v, S){

   # Define this function like this:
   #
   #     p <- nrow(S)
   #     stopifnot(v > 2*p)
   #     v.Anderson <- v - p - 1
   #     return(MCMCpack::riwish(v.Anderson, S))
   #
   # Then these calls are equivalent:
   #
   #     bayessource::riwish_Press(v, S)
   #     W.sources.exact[[i]] <- solve(rWishart(1, v.Anderson, solve(S))[,,1])
   #
   #

   p <- nrow(S)
   stopifnot(v > 2*p)
   v_Anderson <- v - p - 1

   # solve(stats::rWishart(1, v, solve(S))[,,1])
   inv_pd(stats::rWishart(1, v_Anderson, inv_pd(S))[,,1])
}


#' Compute the inverse of a positive-definite matrix
#'
#' Compute the inverse of a positive-definite matrix using Cholesky decomposition.
#'
#' @param X a positive definite matrix
#' @return the inverse of X
#' @family math functions
#' @export
inv_pd <- function(X) {
   chol2inv(chol(X))
}
