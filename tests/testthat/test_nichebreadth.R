library(EcoSimR)
context("Niche overlap null model tests")

tmat <- matrix(rnorm(100),ncol=10,nrow=10)

tmat[c(1,19,13,23,42,90)] <- 0
uval <- length(unique(tmat))



test_that("ra1 algorithm works:",{
  expect_true(is.matrix(ra1()))
  expect_true(is.matrix(ra1(macwarb[,2:5])))
  
}
)


test_that("ra2 algorithm works:",{
  expect_true(is.matrix(ra2()))
  expect_true(is.matrix(ra2(as.matrix(macwarb[,2:5]))))
  expect_equal(length(which(ra2(tmat)==0)),6)
}
)

test_that("ra3 algorithm works:",{
  expect_true(is.matrix(ra3()))
  expect_true(is.matrix(ra3(macwarb[,2:5])))
  expect_equal(length(unique(ra3(tmat))),uval)  
}
)


test_that("ra4 algorithm works:",{
  expect_true(is.matrix(ra4()))
  expect_true(is.matrix(ra4(macwarb[,2:5])))
  expect_equal(sum(apply(ra4(tmat),2,function(x){sum(x==0)}) - apply(tmat,2,function(x){sum(x==0)})), 0)
  expect_equal(sum(apply(ra4(tmat),1,function(x){sum(x==0)}) - apply(tmat,1,function(x){sum(x==0)})), 0)
}
)





