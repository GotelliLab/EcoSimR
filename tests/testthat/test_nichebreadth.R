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

test_that("pianka metric works:"{
expect_true(is.numeric(pianka(ra1())))
expect_true(is.numeric(pianka(ra4(macwarb[,2:5]))))
expect_true(pianka(ra4(macwarb[,2:5])) > .5 )

})

test_that("czekanowski metric works:"{
  expect_true(is.numeric(czekanowski(ra1())))
  expect_true(is.numeric(czekanowski(ra4(macwarb[,2:5]))))
  expect_true(czekanowski(ra4(macwarb[,2:5])) > .5 )
  
})

test_that("pianka_var metric works:"{
  expect_true(is.numeric(pianka_var(ra1())))
  expect_true(is.numeric(pianka_var(ra4(macwarb[,2:5]))))
  expect_true(pianka_var(ra4(macwarb[,2:5])) < .5 )
  
})


test_that("czekanowski_var metric works:"{
  expect_true(is.numeric(czekanowski_var(ra1())))
  expect_true(is.numeric(czekanowski_var(ra4(macwarb[,2:5]))))
  expect_true(czekanowski_var(ra4(macwarb[,2:5])) < .5 )
  
})

test_that("pianka_skew metric works:"{
  expect_true(is.numeric(pianka_skew(ra1())))
  expect_true(is.numeric(pianka_skew(ra4(macwarb[,2:5]))))
  
})


test_that("czekanowski_skew metric works:"{
  expect_true(is.numeric(czekanowski_skew(ra1())))
  expect_true(is.numeric(czekanowski_skew(ra4(macwarb[,2:5]))))  
})







