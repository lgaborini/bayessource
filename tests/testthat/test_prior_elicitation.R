# Test prior elicitation

library(testthat)
library(bayessource)

context('Prior elicitation functions')

# Data ----------------------------------------------------------------------

seed <- 123

# Contains some data.frames with generating results
load('test_dataset_p=2.RData')

# Integer items
# df.background


# Character items
df.background.str <- df.background
df.background.str$source <- as.character(df.background.str$source)

# data.frame and tibbles --------------------------------------------------


test_that('data.frames, integer items work', expect_silent(make_priors_and_init(df.background, col.variables, col.item)))
test_that('data.frames, integer items work', expect_type(make_priors_and_init(df.background, col.variables, col.item), 'list'))

test_that('data.frames, character items work', expect_silent(make_priors_and_init(df.background.str, col.variables, col.item)))
test_that('data.frames, character items work', expect_type(make_priors_and_init(df.background.str, col.variables, col.item), 'list'))


test_that('tibbles, integer items work', {
  skip_if_not_installed('tibble')
  tbl_background <- tibble::as_tibble(df.background)
  expect_silent(make_priors_and_init(tbl_background, col.variables, col.item))
  expect_type(make_priors_and_init(tbl_background, col.variables, col.item), 'list')
})

test_that('tibbles, character items work', {
  skip_if_not_installed('tibble')
  tbl_background.str <- tibble::as_tibble(df.background.str)
  expect_silent(make_priors_and_init(tbl_background.str, col.variables, col.item))
  expect_type(make_priors_and_init(tbl_background.str, col.variables, col.item), 'list')
})



# Parameter checking ------------------------------------------------------


