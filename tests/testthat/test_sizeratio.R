library(EcoSimR)
context("Size ratio null model tests")
test_that("uniform_size algorithm works:",{
  expect_true(is.vector(uniform_size()))
  expect_true(is.vector(uniform_size(rodents$Sonoran)))
  
  }
)

test_that("uniform_size_user algorithm works:",{
  expect_true(is.vector(uniform_size_user()))
  expect_true(is.vector(uniform_size_user(rodents$Sonoran)))
  expect_true(is.vector(do.call(uniform_size_user,list(v = rodents$Sonoran,user.low=4,user.high=20))))
  expect_true(min(do.call(uniform_size_user,list(v = rodents$Sonoran,user.low=4,user.high=20))) > 4)
  expect_true(max(do.call(uniform_size_user,list(v = rodents$Sonoran,user.low=4,user.high=20))) < 20)
  
   
}
)
