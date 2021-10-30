# Use the iris data
head(iris, 3)
col_variables <- c(1,3)
col_item <- 5

# Elicitation using MLE
priors_init <- make_priors_and_init(
   df.background = iris,
   col.variables = col_variables,
   col.item = col_item
)

priors_init

priors_init_2 <- make_priors_and_init(
   df.background = iris,
   col.variables = col_variables,
   col.item = col_item,
   use.priors = "vague",
   alpha_init = 100
)

priors_init_2
