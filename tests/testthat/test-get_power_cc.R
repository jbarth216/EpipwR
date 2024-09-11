test_that("Case Control Output check", {
  expect_s3_class(get_power_cc(100,1000,c(50,100),.05,c(.01,.02,.03)), "data.frame")
})

test_that("Continuous Output check", {
  expect_s3_class(get_power_cc(100,1000,c(50,100),.05,c(.1,.2,.3)), "data.frame")
})


