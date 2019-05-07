#' Elicit priors and initialization from background dataset
#'
#' Elicit priors and initialization from background dataset.
#'
#' Here we fix the hyperparameters for priors on $\theta_i$ and $W_i$, i.e., $B$, $U$, $\mu$ and $n_w$.
#'
#' An appropriate initialization value for $W^{-1}_i$, $i=1,2$ is also generated.
#' Notice that we have three chains in the LR computations, the same initialization is used thrice.
#'
#' ### Priors
#'
#' - use.priors = 'ML': maximum likelihood estimation
#' - use.priors = 'vague': low-information priors
#'    (U is alpha*diag(p), B is beta*diag(p), mu is mu0)
#'
#' The Wishart dofs `nw` are set as small as possible without losing full rank:
#'    $$nw = 2*(p + 1) + 1$$
#'
#' ### Initialization
#'
#' - use.init = 'random': initialize according to the model
#' - use.init = 'vague': low-information initialization (W_i is constant(alpha) + beta*diag(p))
#'
#' Some constants can be changed by passing the new values to ... :
#'
#' - use.priors = 'vague':
#'    alpha = 1, beta = 1, mu0 = 0
#' - use.init = 'vague':
#'    alpha = 1, beta = 100
#'
#' ### Returns
#'
#' A list of variables:
#'
#' - `mu`: the global mean
#' - `B.inv`: the between prior covariance, as the inverse
#' - `U`: the Inverted Wishart scale matrix
#' - `nw`: the Inverted Wishart dof
#' - `W.inv.1`, `W.inv.2`: the within cov initializations, as inverses
#'
#'
#' @param df.background the background dataset
#' @param col.variables columns with variables
#' @param col.item column with item id
#' @param use.priors see details
#' @param use.init see details
#' @param ... additional variables for priors, init
#' @return a list of variables
#' @export
make_priors_and_init <- function(df.background, col.variables, col.item, use.priors = 'ML', use.init = 'random', ...) {

   stopifnot(use.priors %in% c('ML', 'vague'))
   stopifnot(use.init %in% c('random', 'vague'))

   p <- length(col.variables)

   # Capture current dots: optional parameters and preset constants
   dots <- list(...)

   # minimum dof s.t. E[W] is defined, with W ~ Wishart (or Inverse Wishart)
   nw.min <- 2*(p + 1) + 1

   # Perform maximum likelihood estimation on the background dataset
   if (use.priors == 'ML') {
      # On background dataset
      WB <- two.level.multivariate.calculate.UC(df.background, col.variables, col.item)

      mu <- t(WB$all.means)
      B.inv <- solve(WB$B)
      nw <- nw.min

      ## ML prior for W
      ## W ~ IW(nw, U)
      ## We have W0 (sample estimate for W?)
      ## Set U prior s.t.
      ## E[W] = U / (nw - 2*(p + 1)) = W0        # Press

      U <- WB$W * (nw - 2*(p + 1))           # mean of Inverse Wishart (Press)

   }

   # Vague priors
   if (use.priors == 'vague') {

      dots.default <- list(alpha = 1, beta = 1, mu0 = 0)
      # Modified optional parameters
      dots.now <- modifyList(dots.default, dots)

      # Generic priors
      # U is the prior on the Within covariance matrix
      U <- dots.now$alpha * diag(p)
      mu <- rep(dots.now$mu0, p)
      # B.inv is the prior Between precision matrix
      B.inv <- (1 / dots.now$beta) * diag(p)

      nw <- nw.min
   }


   # End of priors
   priors <- list()
   priors$mu <- mu
   priors$U <- U
   priors$B.inv <- B.inv
   priors$nw <- nw

   #' ## Chain initialization
   #'
   #' Here we initialize the chain by supplying $W_1^{-1}$ and $W_2^{-1}$.
   #' The $\theta_i$ are generated accordingly to the model.
   #'
   #' Three methods:
   #'
   #' - sample from the hierarchical model using the hyperparameters $n_w$, $U$
   #' - initialize from fixed values
   #' - cheat and use the known values

   # Initialize according to the model
   # Warning: RNG is used here!
   if (use.init == 'random') {

      # W.1 <- riwish(nw.exact, U)
      W.1 <- riwish_Press(nw, U)
      W.inv.1 <- solve(W.1)

      # is there a better option?
      W.inv.2 <- W.inv.1
      rm(W.1)
   }

   # Initialize at a fixed value
   if (use.init == 'vague') {

      dots.default <- list(alpha = 1, beta = 100)
      dots.now <- modifyList(dots.default, dots)

      # W.1 is a covariance matrix, not a precision!
      W.1 <- dots.now$alpha*matrix(1, p, p) + dots.now$beta*diag(p)
      W.inv.1 <- solve(W.1)
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
#' Returns minimum value for degrees of freedom such that Inverted Wishart has a mean..
#'
#' Uses Press[1] parametrization.
#'
#' \deqn{X ~ IW(nw, S)}, with \eqn{S = pxp} matrix, \eqn{nw > 2p} (the degrees of freedom).
#' Then:
#' \deqn{E[X] = S / (nw - 2(p + 1))}
#'
#' Finally, the minimum dof $v$ is $ nw = 2*(p + 1) + 1 $.
#'
#' @references [1]S. J. Press, Applied multivariate analysis: using Bayesian and frequentist methods of inference. Courier Corporation, 2012.
#'
#' @param p dimension
#' @return minimum dof
#' @export
get_minimum_nw_IW <- function(p) {
   if (p < 1) { stop('p must be >= 1')}

   nw.min <- 2*(p + 1) + 1
   nw.min
}



#' Generate random sample from Inverted Wishart.
#'
#' Using Press parametrization.
#'
#' Uses \pkg{MCMCpack}::riwish.
#'
#' @param v dof (\eqn{> 2p})
#' @param S the scale matrix (pxp)
#' @return a single random variate from IW(v, S)
#' @template InverseWishart_Press
riwish_Press <- function(v, S){
   p <- nrow(S)
   stopifnot(v > 2*p)
   v.Anderson <- v - p - 1
   return(MCMCpack::riwish(v.Anderson, S))
}
