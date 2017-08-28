context("Test function initializeGphmm.")

test_that("Function initializeGphmm", {
  param = initializeGphmm()
  expect_true(length(param) == 10)
  expect_true(class(param$pp) == 'matrix')
  expect_true(sum(colSums(param$pp)) == 4)
  expect_true(sum(rowSums(param$pp)) == 4)
  expect_true(sum(param$qR) == 1)
  expect_true(sum(param$qX) == 1)
  expect_true(sum(param$qY) == 1)
  expect_true(length(param$deltaX) == 2)
  expect_true(length(param$deltaY) == 2)
})
