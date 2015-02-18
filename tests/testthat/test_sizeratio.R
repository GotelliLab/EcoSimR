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
  expect_true(is.vector(do.call(uniform_size_user,list(speciesData = rodents$Sonoran,userLow=4,userHigh=20))))
  expect_true(min(do.call(uniform_size_user,list(speciesData = rodents$Sonoran,userLow=4,userHigh=20))) > 4)
  expect_true(max(do.call(uniform_size_user,list(speciesData = rodents$Sonoran,userLow=4,userHigh=20))) < 20)
  }
)

test_that("source_pool_draw algorithm works:",{
  expect_true(is.vector(source_pool_draw()))
  expect_true(is.vector(source_pool_draw(rodents$Sonoran)))
  expect_true(is.vector(do.call(source_pool_draw,list(speciesData = rodents$Sonoran,
       sourcePool = runif(1000,min(rodents$Sonoran),max(rodents$Sonoran)),
       speciesProbs = rbeta(1000,1,1) ))))
}
)

test_that("gamma_size algorithm works:",{
  expect_true(is.vector(size_gamma()))
  expect_true(is.vector(size_gamma(rodents$Sonoran))) 
}
)


test_that("min_diff metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(min_diff()))
  ### Test that minimum difference is accurate
  expect_equal(min_diff(c(1,1.2,-3,3)),.2) 
})

test_that("min_diff metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(min_diff()))
  ### Test that minimum difference is accurate
  expect_equal(min_diff(c(1,1.2,-3,3)),.2) 
  expect_equal(min_diff(c(1,1.2,1.1,-3,3)),.1) 
  
  
  })

test_that("min_ratio metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(min_ratio()))
  ### Test that minimum ratio is accurate
  expect_equal(min_ratio(c(1,2,3,4,5,6)),(6/5)) 
  
})


test_that("var_diff metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(var_diff()))
})


test_that("var_ratio metric works",{
  ### Test that proper object is returned
  expect_true(is.numeric(var_ratio()))
  })

test_that("size_null_model works with all combinations of metrics and algorithms",{
  ### Test that proper object is returned
  expect_is(size_null_model(rodents,metric ="Min.Diff" ,algo = "Uniform.Size",nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Min.Ratio" ,algo = "Uniform.Size",nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Diff" ,algo = "Uniform.Size",nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Ratio" ,algo = "Uniform.Size",nRep=10),"sizenullmod")
  
  expect_is(size_null_model(rodents,metric ="Min.Diff" ,algo = "Uniform.Size.User",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Min.Ratio" ,algo = "Uniform.Size.User",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Diff" ,algo = "Uniform.Size.User",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Ratio" ,algo = "Uniform.Size.User",algoOpts = list(userLow = 3,userHigh=15),nRep=10),"sizenullmod")
  
  
  expect_is(size_null_model(rodents,metric ="Min.Diff" ,algo = "Source.Pool",algoOpts = list(sourcePool = runif(1000,min(rodents$Sonoran),max(rodents$Sonoran)),
                                                                                             speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Min.Ratio" ,algo = "Source.Pool",algoOpts = list(sourcePool = runif(1000,min(rodents$Sonoran),max(rodents$Sonoran)),
                                                                                              speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Diff" ,algo = "Source.Pool",algoOpts = list(sourcePool = runif(1000,min(rodents$Sonoran),max(rodents$Sonoran)),
                                                                                                      speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Ratio" ,algo = "Source.Pool",algoOpts = list(sourcePool = runif(1000,min(rodents$Sonoran),max(rodents$Sonoran)),
                                                                                                       speciesProbs = rbeta(1000,1,1) ),nRep=10),"sizenullmod")
   
  expect_is(size_null_model(rodents,metric ="Min.Diff" ,algo = "Gamma",nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Min.Ratio" ,algo = "Gamma",nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Diff" ,algo = "Gamma",nRep=10),"sizenullmod")
  expect_is(size_null_model(rodents,metric ="Var.Ratio" ,algo = "Gamma",nRep=10),"sizenullmod")
  
  smod <- size_null_model(rodents,metric ="Var.Diff" ,algo = "Gamma",nRep=10)
  
  expect_output(summary(smod),"Metric:  Var.Diff")
  
})

