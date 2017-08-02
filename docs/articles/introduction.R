## ---- include = FALSE----------------------------------------------------
library(bayessource)
knitr::opts_chunk$set(fig.datapi = 96)

## ---- results = FALSE----------------------------------------------------
set.seed(123)

p <- 2
mean.quest <- c(0, 0)
mean.ref <- c(3, 2)
cov.quest <- rwish(3, diag(2))
cov.ref <- rwish(5, diag(2))
n.quest <- 20
n.ref <- n.quest

df.quest <- data.frame(rmvnorm(n.quest, mean.quest, cov.quest))
df.ref <- data.frame(rmvnorm(n.ref, mean.ref, cov.ref))

## ------------------------------------------------------------------------

library(ggplot2)
ggplot() + 
   geom_point(aes(x = X1, y = X2), col = 'red', data = df.quest) +
   geom_point(aes(x = X1, y = X2), col = 'blue', data = df.ref)

## ------------------------------------------------------------------------
eps <- 0.001
B.inv <- eps*diag(p)
W.inv.1 <- eps*diag(p)
W.inv.2 <- eps*diag(p)
U <- eps*diag(p)
nw <- 2*(p + 1) + 1
mu <- (mean.quest + mean.ref)/2

## ----collapse=TRUE-------------------------------------------------------
burn.in = 1000
n.iter = 10000
marginalLikelihood(as.matrix(df.quest), n.iter, B.inv, W.inv.1, U, nw, mu, burn.in, verbose = FALSE)

## ----collapse=TRUE-------------------------------------------------------
samesource_C(as.matrix(df.quest), as.matrix(df.ref), n.iter, B.inv, W.inv.1, W.inv.2, U, nw, mu, burn.in, verbose = FALSE)

## ----collapse=TRUE-------------------------------------------------------
samesource_C(as.matrix(df.ref)[1:20,], as.matrix(df.ref), n.iter, B.inv, W.inv.1, W.inv.2, U, nw, mu, burn.in, verbose = FALSE)

## ------------------------------------------------------------------------
results <- marginalLikelihood(as.matrix(df.quest), n.iter, B.inv, W.inv.1, U, nw, mu, burn.in, output.mcmc = TRUE)

## ----collapse=TRUE-------------------------------------------------------
head(results$mcmc, 4)

## ------------------------------------------------------------------------
library(coda)
summary(results$mcmc)

## ----fig.height=3--------------------------------------------------------
traceplot(results$mcmc)

## ---- collapse=TRUE------------------------------------------------------
n.samples <- nrow(results$mcmc)
W.inv.samples <- results$mcmc[, paste0('W.inv.', seq(1:(p^2)))]
head(W.inv.samples, 5)
W.inv.samples.cube <- array(W.inv.samples, dim = c(n.samples, p, p))
dim(W.inv.samples.cube)

