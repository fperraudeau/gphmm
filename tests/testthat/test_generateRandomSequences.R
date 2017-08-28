context("Test function generateRandomSequences.")

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
