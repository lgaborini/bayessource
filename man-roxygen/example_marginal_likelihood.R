# Use the iris data
df_X <- iris[, c(1,3)]
X <- as.matrix(df_X)
p <- ncol(X)

# Observations
plot(X)
# True species
plot(X, col = iris$Species)


# Priors
B.inv <- diag(p)
W.inv <- diag(p)
U <- diag(p)
mu <- colMeans(X)
# d.o.f.
nw <- get_minimum_nw_IW(p)

# Computational parameters
n_iter <- 10000
n_burn_in <- 1000

# Compute the marginal likelihood
marginalLikelihood(
   X = X,
   n.iter = n_iter,
   B.inv = B.inv,
   W.inv = W.inv,
   U = U,
   nw = nw,
   burn.in = n_burn_in,
   mu = mu
)

# Diagnostics ------------

# Retain the full Gibbs output
list_mcmc <- marginalLikelihood(
   X = X,
   n.iter = 10000,
   B.inv = B.inv,
   W.inv = W.inv,
   U = U,
   nw = get_minimum_nw_IW(p),
   burn.in = 1000,
   mu = mu,
   output.mcmc = TRUE
)

# The log-ML value
list_mcmc$value

# The full Gibbs chain output
library(coda)
head(list_mcmc$mcmc, 20)

# Diagnostics plots by coda

plot(list_mcmc$mcmc)
coda::autocorr.plot(list_mcmc$mcmc)
