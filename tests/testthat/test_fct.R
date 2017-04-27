context("Test functions individually.")

test_that("Function generateRandomSequences", {
  seq = generateRandomSequences(1, 5, 0, 5373)
  expect_true(as.character(seq) == 'TAGAA')
  expect_true(width(seq) == 5)
  expect_true(length(seq) == 1)
  expect_true(class(seq) == 'DNAStringSet')
  
  seqs = generateRandomSequences(2, 5, 0, 5373)
  expect_true(as.character(seqs)[1] == 'TAGAA')
  expect_true(width(seqs)[1] == 5)
  expect_true(length(seqs) == 2)
  expect_true(class(seqs) == 'DNAStringSet')
})


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


test_that("Function initializeGphmm", {
  param = initializeGphmm()
  prob = computegphmm('ATG', 'ATGC', param, 20, 'short')
  expect_that(round(prob), equals(-39))
  prob = computegphmm('ATG', 'ATGC', param, 20, 'long')
  expect_that(prob$path, equals('MMMD'))
  expect_that(round(prob$V), equals(-39))
  expect_that(length(prob), equals(4))
  expect_that(names(prob), equals(c('V', 'path', 'read', 'ref')))
})
