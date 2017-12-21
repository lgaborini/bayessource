# Run test case
#
# Actually check that code does not throw
#

rm(list = ls())

library(bayessource)
library(testthat)
context('samesource_C')

seed <- 123

# Generate data ------------------------

load('test_dataset_p=2.RData')

# Make bayessource-friendly datasets
ref <- df.ref[, col.variables] %>% as.matrix
quest.same <- df.quest.same[, col.variables] %>% as.matrix
quest.diff <- df.quest.diff[, col.variables] %>% as.matrix
quest.ref.same <- rbind(quest.same, ref)
quest.ref.diff <- rbind(quest.diff, ref)

# Use exact parameters
theta.H1.exact <- list.exact$theta.sources.exact[[1]]
W.H1.exact <- list.exact$W.sources.exact[[1]]
theta.H2.exact <- do.call(rbind, list.exact$theta.sources.exact[list.exact$source.quest])
W.H2.exact <- list.exact$W.sources.exact[[list.exact$source.quest]]

# Set-up short simulation
burn.in <- 100
n.iter <- 1000

B.inv <- solve(list.exact$B.exact)
W.inv.1 <- solve(W.H1.exact)
W.inv.2 <- solve(W.H2.exact)
U <- list.exact$U.exact
nw <- list.exact$nw.exact
mu <- list.exact$mu.exact


# H1 ----------------------------------------------------------------------

test_that('C code does not throw', testthat::expect_silent(samesource_C(quest.same, ref, 100, B.inv, W.inv.1, W.inv.2, U, nw, mu, 10)))

set.seed(seed)
f.LLR.same <- samesource_C(quest.same, ref, n.iter, B.inv, W.inv.1, W.inv.2, U, nw, mu, burn.in)
f.LLR.same

test_that('C code returns LR > 1 under H1 (same sources)', expect_gte(f.LLR.same, 0))

# H2 ----------------------------------------------------------------------

set.seed(seed)
f.LLR.diff <- samesource_C(quest.diff, ref, n.iter, B.inv, W.inv.1, W.inv.2, U, nw, mu, burn.in)
f.LLR.diff

test_that('C code returns LR < 1 under H2 (different sources)', expect_gte(-f.LLR.diff, 0))

