library(testthat)
library(igraph)

test_that("autoClusterMapper errors if no predicted$class", {
  obj <- methods::new("TDAobj")
  obj@graph <- list(igraph = make_ring(3), predicted = NULL)
  expect_error(
    autoClusterMapper(obj),
    "`x@graph\\$predicted\\$class` is not set"
  )
})

test_that("autoClusterMapper errors if igraph missing", {
  obj <- methods::new("TDAobj")
  obj@graph <- list(igraph = NULL, predicted = list(class="ER"))
  expect_error(autoClusterMapper(obj),
               "`x@graph\\$igraph` is not set")
})
