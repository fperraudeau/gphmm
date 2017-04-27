Sys.setenv("R_TESTS" = "")
library(testthat)
library(gphmm)

test_check("gphmm")
