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
#' Multivariate functions supporting paper:
#' Bozza, A., Taroni, F., Marquis, R., M. Schmittbuhl
#'
#' This function is RNG-free.
#'
#' @param population matrix or dataframe containing data columns, item column and other columns
#' @param variables data columns (indexes or character vectors)
#' @param idx.item item column (index or character vector)
#'
#' @return list with items:
#' - group.means (matrix, row*variable)
#' - all.means (vector)
#' - W (the within covariance matrix)
#' - B (the between covariance matrix)
#'
two.level.multivariate.calculate.UC <- function(population, idx.variables, idx.item) {

   # define how many unique items are in the population
   all.items <- unique(population[, idx.item])
   n.items <- length(all.items)

   variable.names <- colnames(population[, idx.variables])
   n.variables <- length(idx.variables)
   if (n.variables < 2) {
      stop("cannot handle fewer than two variables with this function")
   }

   # How many replicates per item
   # table(population[, idx.item])

   # define the output matrices
   Sw <- matrix(0, ncol = n.variables, nrow = n.variables)
   rownames(Sw) <- variable.names
   colnames(Sw) <- variable.names
   S <- Sw

   all.means <- matrix(apply(population[, idx.variables], 2, mean), nrow = 1)
   colnames(all.means) <- variable.names

   # dummy within group means
   # rows: items
   # columns: variables
   group.means <- matrix(0, nrow = n.items, ncol = n.variables)
   rownames(group.means) <- as.character(all.items)
   colnames(group.means) <- variable.names

   n.pop <- dim(population)[1]
   # for each unique item
   for (i.item in 1:n.items) {

      # which rows refer to repeated measurements on the same item
      idx.replicates <- which(population[, idx.item] == all.items[i.item])
      n.replicates <- length(idx.replicates)

      # pick out the measurements for the item
      tempdat <- as.matrix(population[idx.replicates, idx.variables])

      # assign values to the group means for the population
      means <- colMeans(tempdat)
      group.means[i.item, ] <- means

      # sum the within and between item means
      S.this <- t(means - all.means) %*% (means - all.means) * n.replicates
      S <- S + S.this

      # sum for all measurements
      Cov.this <- cov(tempdat)*(n.replicates - 1)
      Sw <- Sw + Cov.this
   }

   # convert matrix Sw to matrix U
   W <- Sw/(n.pop - n.items)

   # convert matrix S to matrix C
   B <- S/(n.items - 1)

   op <- list(group.means = group.means, all.means = all.means, W = W, B = B)
   return(op)
}

