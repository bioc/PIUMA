library(testthat)

test_that("scaleData_01 rescales correctly", {
  expect_equal(scaleData_01(c(2,4)), c(0,1))
})

test_that("jaccard gives correct ratio", {
  expect_equal(jaccard(1:3, 2:4), 0.5)
})

test_that("Plotting function existence", {
     expect_true(is.function(plot_ScaleFreeLaw))
     expect_true(is.function(plot_projection_plot))
})

