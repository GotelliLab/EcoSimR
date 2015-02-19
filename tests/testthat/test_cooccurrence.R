library(EcoSimR)
context("Co-occurrence null model tests")
tmat <- matrix(rbinom(100,1,.4),ncol=10,nrow=10)
realData <- as.matrix(wiFinches[,2:20])
testthat("vector_sample algorithm works",{
  expect_true(is.vector(vector_sample(speciesData=rbinom(10,1,0.5),weights=runif(10))))
})

testthat("sim1 algorithm works",{
  expect_true(is.matrix(sim1(tmat)))
  expect_true(is.matrix(sim1(realData)))
})

testthat("sim2 algorithm works",{
  expect_true(is.matrix(sim2(tmat)))
  expect_true(is.matrix(sim2(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim2(tmat),1,sum),apply(tmat,1,sum))
})

testthat("sim3 algorithm works",{
  expect_true(is.matrix(sim3(tmat)))
  expect_true(is.matrix(sim3(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim3(tmat),2,sum),apply(tmat,2,sum))
})


testthat("sim4 algorithm works",{
  expect_true(is.matrix(sim4(tmat)))
  expect_true(is.matrix(sim4(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim4(tmat),1,sum),apply(tmat,1,sum))
})

testthat("sim5 algorithm works",{
  expect_true(is.matrix(sim5(tmat)))
  expect_true(is.matrix(sim5(realData)))
  ### Test that row sums are preserved
  expect_equal(apply(sim5(tmat),2,sum),apply(tmat,2,sum))
})
  
testthat("sim6 algorithm works",{
  expect_true(is.matrix(sim6(tmat)))
  expect_true(is.matrix(sim6(realData)))
})



