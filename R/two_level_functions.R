# This script requires no external library.
#
#
# Calculation of ML Between/Within covariance matrices for a two-level hierarchical model.
#
#------------------------------------------------------------

#' Calculate Between/Within covariance matrices
#'
#' Two-level function: item is chosen with the third parameter.
#'
#' Items who are observed only once cannot contribute to within covariance estimation.
#' In that case, either we skip its contribute (default) or fail.
#'
#' Multivariate functions supporting paper:
#'
#' Bozza, A., Taroni, F., Marquis, R., M. Schmittbuhl - "Probabilistic evaluation of handwriting evidence: likelihood ratio for authorship"
#'
#' @param population matrix or dataframe containing data columns, item column and other columns
#' @param idx.variables data columns (indexes or character vectors)
#' @param idx.item item column (index or character vector)
#' @param singular.action action to take when one of the items is observed once; one of `'skip'` (default) or `'fail'`
#' @importFrom stats cov
#' @importFrom utils modifyList
#' @return list with items:
#' - `group.means` (matrix, row*variable)
#' - `all.means` (vector)
#' - `W` (the within covariance matrix)
#' - `B` (the between covariance matrix)
#' @keywords internal
#' @family core functions
two.level.multivariate.calculate.UC <- function(population, idx.variables, idx.item, singular.action = 'skip') {

   stopifnot(singular.action %in% c('skip', 'fail'))

   # Convert to matrix
   mtx_population <- as.matrix(population[, idx.variables])

   # Convert items to factor
   item_population <- as.matrix(population[, idx.item])

   # define how many unique items are in the population
   all.items <- as.vector(unique(item_population))
   n.items <- length(all.items)

   if (n.items == 1) {
      stop('two.level.multivariate.calculate.UC: only one item in population')
   }

   variable.names <- colnames(population)[idx.variables]
   n.variables <- length(idx.variables)
   if (n.variables < 2) {
      stop("two.level.multivariate.calculate.UC: cannot handle fewer than two variables with this function")
   }

   # How many replicates per item
   # table(population[, idx.item])

   # define the output matrices
   Sw <- matrix(0, ncol = n.variables, nrow = n.variables)
   rownames(Sw) <- variable.names
   colnames(Sw) <- variable.names
   S <- Sw

   all.means <- matrix(colMeans(mtx_population), nrow = 1)
   # all.means <- matrix(apply(population[, idx.variables], 2, mean), nrow = 1)
   colnames(all.means) <- variable.names

   # dummy within group means
   # rows: items
   # columns: variables
   group.means <- matrix(0, nrow = n.items, ncol = n.variables)
   rownames(group.means) <- as.character(all.items)
   colnames(group.means) <- variable.names

   n.pop <- dim(population)[1]
   n.items.skip <- 0

   # for each unique item
   for (i.item in 1:n.items) {

      this_item <- all.items[i.item]

      # which rows refer to repeated measurements on the same item
      idx.replicates <- which(item_population == this_item)
      n.replicates <- length(idx.replicates)

      # pick out the measurements for the item
      tempdat <- mtx_population[idx.replicates, , drop = FALSE]


      # assign values to the group means for the population
      means <- colMeans(tempdat)
      group.means[i.item, ] <- means

      if (n.replicates == 1) {
         if (identical(singular.action, 'skip')) {
            message(paste0('two.level.multivariate.calculate.UC: item "', i.item, '" has 1 replicate: cannot contribute to within covariance matrix estimation. Skipping.'))
            n.items.skip <- n.items.skip + 1
            next
         }
         stop(paste0('two.level.multivariate.calculate.UC: item "', i.item, '" has 1 replicate: cannot contribute to within covariance matrix estimation. Failing.'))
      }

      # sum the within and between item means
      S.this <- t(means - all.means) %*% (means - all.means)
      S <- S + S.this

      # Covariance for this item, without normalization factor
      # = sum of squares to the item means
      #
      # sum for all measurements
      Cov.this <- cov(tempdat)*(n.replicates - 1)
      Sw <- Sw + Cov.this
   }

   # Only these items have contributed to S, Sw
   n.items.contributing <- n.items - n.items.skip

   # convert matrix Sw to matrix U
   W <- Sw/(n.pop - n.items.contributing)

   # convert matrix S to matrix C
   B <- S/(n.items.contributing - 1)

   op <- list(group.means = group.means, all.means = all.means, W = W, B = B)
   return(op)
}

