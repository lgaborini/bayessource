#
#
# Improved R code
#

# source('statistical_functions.R')

#' Bayesian same source hypothesis. Gaussian MV.
#' Implemented in R (slow).
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
#' @return a LR value, or a list(LR value, posterior samples)
#' @export
#' @template gaussmv_model
samesource_R <- function(dati, n.iter, B.inv, W.inv, U, nw, mu, burn.in, output.mcmc = FALSE, verbose = FALSE) {

   # Redundant arguments
   p <- nrow(B.inv)
   # vp <- 1:p
   # lc0 <- -(nw - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((nw - p - vp)/2)))

   # \theta chain outputs
   theta.gibbs <- list()
   B.upd.gibbs <- list()
   mu.upd.gibbs <- list()

   # W chain outputs
   W.inv.gibbs <- list()
   # W.gibbs <- list()
   U.upd.gibbs <- list()

   nr <- dim(dati)[1]

   bary <- matrix(colMeans(dati), nrow = 1)

   ## Original
   # S <- 0
   # for (j in 1:nr) {
   #    S <- S + t(dati[j, ] - bary) %*% (dati[j, ] - bary)
   # }

   ## Centered
   # dati.centered <- sweep(dati, 2, bary, '-')
   dati.centered <- scale(dati, center = bary, scale = FALSE)
   S <- crossprod(dati.centered)

   ## By definition:
   # S <- (nrow(dati) - 1) * cov(dati)

   nw.star <- nr + nw
   B.inv.mu <- B.inv %*% mu

   # Initialize the Gibbs sampler
   # Gibbs chain outputs: W.inv.g, theta.g
   W.inv.g <- W.inv
   # B.inv.g <- B.inv
   # mu.g <- mu
   # U.g <- U

   # Maximum likelihood
   logf.star <- -Inf
   detW.star <- NA

   for (i in 1:n.iter) {
      ### Sample theta:
      # \theta ~ N(\mu, B)
      # Update \theta prior covariance B
      # B.upd.inv.g <- B.inv + nr * W.inv.g
      B.upd.g <- solve(B.inv + nr * W.inv.g)
      # Update \theta prior mean \mu
      mu.upd.g <- B.upd.g %*% (B.inv.mu + nr * W.inv.g %*% t(bary))
      # Sample \theta from N(\mu, B)
      #theta <- t(mvrnorm(1, mu.upd.g, B.upd.g))  # as row vector
      theta.g <- matrix(MASS::mvrnorm(1, mu.upd.g, B.upd.g), nrow = 1)

      # Save chain outputs
      theta.gibbs[[i]] <- theta.g
      B.upd.gibbs[[i]] <- B.upd.g
      mu.upd.gibbs[[i]] <- mu.upd.g

      ### Sample W and its inverse
      # W ~ IW(U, nw)
      # Update W's scale matrix U
      # U.upd.g <- nr * t(theta.g - bary) %*% (theta.g - bary) + U + S
      U.upd.g <- nr * crossprod(theta.g - bary) + U + S
      # Note: rank-1 update to U+S
      # det(U.upd.g) can be easily computed (to check!)
      # det(U.upd.g) = det(U + S)*(1 + nr^2 * (theta.g - bary) %*% t(theta.g - bary))

      ## Explicit Inverse Wishart sampling
      # W.g <- riwish(nw.star, U.upd.g)
      # W.inv.g <- solve(W.g)
      # detW <- det(W.g)

      ## Using Wishart definition: avoid generating (and inverting) W
      ## riwish(v, S) = solve(rwish(v,solve(S)))
      # W.inv.g <- rwish(nw.star, solve(U.upd.g))

      ## Using base::R C function instead of R MCMCpack (faster)
      W.inv.g <- rWishart(1, nw.star, solve(U.upd.g))[,,1]
      detW.g <- 1/det(W.inv.g)

      # Save chain outputs
      U.upd.gibbs[[i]] <- U.upd.g
      W.inv.gibbs[[i]] <- W.inv.g
      # W.gibbs[[i]] <- solve(W.inv.g)

      # Update posteriors -> priors
      # mu.g <- mu.upd.g
      # U.g <- U.upd.g
      # B.inv.g <- solve(B.upd.g)
      # B.inv.g <- B.upd.inv.g

      # dall'output del Gibbs sampling posso trovarmi psi_star
      # calcolo la densità delle mie osservazioni in ogni punto
      # della catena e prendo la moda
      # theta sarà quello dell'iterazione corrente
      # W sarà quello dell'iterazione corrente

      # logf0 <- -(p * nr/2) * log(2 * pi) - (nr/2) * log(detW.g)

      # dati.centered.theta <- sweep(dati, 2, theta.g, '-')
      #
      # logf <- logf0
      # for (j in 1:nr) {
      #    logf <- logf - (dati.centered.theta[j,] %*% (W.inv.g %*% dati.centered.theta[j,]))/2
      # }
      #
      # logf <- logf0
      # for (j in 1:nr) {
      #    logf <- logf - (dati[j, ] - theta.g) %*% W.inv.g %*% t(dati[j, ] - theta.g)/2
      # }

      # Faster:
      # logf <- logf0 - sum(diag(dati.centered.theta %*% (W.inv.g %*% t(dati.centered.theta))))/2


      # By definition:
      logf <- sum(dmvnorm_fast(dati, theta.g, solve(W.inv.g), log = TRUE, is.chol = FALSE))

      if (logf > logf.star) {
         if (verbose) {
            print(sprintf('New minimum found: logf = %g, previous = %g, delta = %g', logf, logf.star, logf - logf.star))
         }

         theta.star <- theta.g
         W.inv.star <- W.inv.g
         logf.star <- logf
         detW.star <- detW.g

      }
   }

   # ML parameter estimates: psi.star
   # {theta.star, W.inv.star}

   rm(B.upd.g, mu.upd.g, theta.g, U.upd.g)

   # save(list = ls(), file = 'samesource_num_end_ML.RData')

###########################################
##### ESTIMATE THE POSTERIOR ORDINATES ####
###########################################

   # Inverse Wishart normalization constant
   # lc0.star <- -(nw.star - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((nw.star - p - (1:p))/2)))

   idx.iter <- (burn.in + 1):n.iter
   lpihat.theta.star.samples <- rep(NA, n.iter)
   lpihat.W.star.samples <- rep(NA, n.iter)
   for (i in idx.iter) {
      # Gibbs samples
      # theta.g <- theta.gibbs[[i]]
      mu.upd.g <- mu.upd.gibbs[[i]]
      B.upd.g <- B.upd.gibbs[[i]]

      # W.inv.g <- W.inv.gibbs[[i]]
      U.upd.g <- U.upd.gibbs[[i]]

      ### \theta block (using W.inv.g => mu.upd.g, B.upd.g)
      # mu.upd.g: la media della condizionata completa di theta
      # B.upd.g: la matrice varianza-cov della condizionata completa di theta

      # B.upd.g <- solve(B.inv + nr * W.inv.g)
      # mu.upd.g <- B.upd.g %*% (invBmu + nr * W.inv.g %*% t(bary))

      # Estimate \theta* posterior density in \theta.star using all W.inv.g
      # lpihat.theta.star <- dmvnorm(t(theta.star), t(mu.upd.g), B.upd.g, log = TRUE)
      lpihat.theta.star <- dmvnorm_fast(theta.star, t(mu.upd.g), B.upd.g, log = TRUE, is.chol = FALSE)
      lpihat.theta.star.samples[i] <- lpihat.theta.star

      ### W block (using theta.g => U.upd.g)
      # U.upd.g: la matrice di scala della condizionata completa di W
      # nw.star: gradi di libertà della condizionata completa di W

      # U.upd.g <- nr * t(theta.g - bary) %*% (theta.g - bary) + U + S
      # U.upd.g <- nr * crossprod(theta.g - bary) + U + S

      # Estimate W* posterior density in W.star using all theta.g
      # Notice that W.star always appears as its inverse:

      # Distribution of inverse:
      # if W.star ~ IWishart(n, U), solve(W.star) = W.inv.star ~ Wishart(n, solve(U))
      # !!! must correct the pdf with the Jacobian !!!

      # lpihat.W.star <- diwishart(solve(W.inv.star), nw.star, U.upd.g, log = TRUE, is.chol = FALSE)          # OK
      # lpihat.W.star <- dwishart(W.inv.star, nw.star, solve(U.upd.g), log = TRUE, is.chol = FALSE)           # !!! must correct with Jacobian !!!

      lpihat.W.star <- diwishart_inverse(W.inv.star, nw.star, U.upd.g, log = TRUE, is.chol = FALSE)          # OK

      # Explicit Inverse Wishart pdf
      # should be equal to diwishart(solve(W.inv.star), nw.star, U.upd.g)
      # lpihat.W.star <- lc0.star + (nw.star - p - 1)/2 * (log(det(U.upd.g))) - sum(diag(W.inv.star %*% U.upd.g))/2 - nw.star/2 * log(detW.star)

      lpihat.W.star.samples[i] <- lpihat.W.star



   }

   # The posterior marginals as MC mean:
   # lpihat.theta.star <- lpihat.theta.star/length(idx.iter)
   # lpihat.W.star <- lpihat.W.star/length(idx.iter)

   # MC mean:
   lpihat.theta.star <- logSumExpMean(lpihat.theta.star.samples[idx.iter])
   lpihat.W.star <- logSumExpMean(lpihat.W.star.samples[idx.iter])

   # W* posterior density: constants + MC mean
   # lDetW.inv.star <- log(det(W.inv.star))
   # lc0star <- -(nw.star - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((nw.star - p - vp)/2)))
   # lpihat.W.star <- lc0star + nw.star * lDetW.inv.star/2 + lpihat.W.star

   # Joint posterior on \theta, W: independent parameters
   lpihat.psi.star <- lpihat.theta.star + lpihat.W.star

   if (verbose) {
      print('Marginal posterior ordinates')
      print(sprintf('  lpihat.theta.star: %g', lpihat.theta.star))
      print(sprintf('  lpihat.W.star: %g', lpihat.W.star))
      # print(sprintf('  lc0star: %g', lc0star))
   }

##########################################
##### COMPUTE THE PRIOR ORDINATES ########
##########################################

   # lpi.theta.star <- dmvnorm(t(theta.star), t(mu), solve(invB), log = TRUE)
   lpi.theta.star <- dmvnorm_fast(theta.star, t(mu), solve(B.inv), log = TRUE, is.chol = FALSE)

   # Explicit Inverse Wishart pdf
   # lc0 <- -(nw - p - 1) * p/2 * log(2) - p * (p - 1)/4 * log(pi) - sum(log(gamma((nw - p - (1:p))/2)))
   # lpi.W.star <- lc0 + (nw - p - 1)/2 * (log(det(U)))  - nw/2 * log(detW.star) - sum(diag(W.inv.star %*% U))/2
   # Inverse Wishart from inverse
   lpi.W.star <- diwishart_inverse(W.inv.star, nw, U, log = TRUE, is.chol = FALSE)

   # lpi.W.star <- diwishart(solve(W.inv.star), nw, U, log = TRUE, is.chol = FALSE)          # !!! must correct with Jacobian !!!
   lpi.psi.star <- lpi.theta.star + lpi.W.star

   if (verbose) {
      print('Marginal prior ordinates')
      print(sprintf('  lpi.theta.star: %g', lpi.theta.star))
      print(sprintf('  lpi.W.star: %g', lpi.W.star))
   }

##########################################
##### OBTAIN THE MARGINAL DENSITY ########
##########################################
   lmhatHp.num <- logf.star + lpi.psi.star - lpihat.psi.star      # Equation (8)

   if (verbose) {
      print('Marginal density')
      print(sprintf('logf.star: %g', logf.star))
      print(sprintf('lpi.psi.star: %g', lpi.psi.star))
      print(sprintf('lpihat.psi.star: %g', lpihat.psi.star))
      print(sprintf('lmhatHp.num: %g', lmhatHp.num))
   }

   # Make and output MCMC object for diagnostics
   if (output.mcmc == TRUE) {
      theta.mtx <- t(sapply(theta.gibbs, cbind))
      colnames(theta.mtx) <- paste0('theta.', 1:p)

      W.inv.mtx <- t(sapply(W.inv.gibbs, as.numeric))
      colnames(W.inv.mtx) <- paste0('W.inv.', 1:(p^2))

      mcmc.data <- cbind(
         theta.mtx,
         W.inv.mtx)

      return(
         list(
            # mcmc = mcmc(data = mcmc.data, start = 1, end = i),
            mcmc = mcmc(data = mcmc.data[(burn.in + 1):n.iter, ], start = (burn.in + 1), end = i),
            LR.num = lmhatHp.num
         )
      )
   }

   # save(list = ls(), file = 'samesource_num.RData')
   return(lmhatHp.num)

}
