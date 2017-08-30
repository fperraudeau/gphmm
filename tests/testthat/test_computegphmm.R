context("Test function computegphmm.")

test_that("Function computegphmm", {
  param = initializeGphmm()
  prob = computegphmm('ATG', 'ATGC', param, qv = 20, output = 'short')
  expect_that(round(prob), equals(-13))
  prob = computegphmm('ATG', 'ATGC', param, 20)
  expect_that(prob$path, equals('MMMD'))
  expect_that(round(prob$V), equals(-13))
  expect_that(length(prob), equals(4))
  expect_that(names(prob), equals(c('V', 'path', 'read', 'ref')))
})
