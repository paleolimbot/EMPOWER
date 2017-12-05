context("model_params")

test_that("zero-length model params are allowed", {
  zlmp <- model_params()
  expect_is(zlmp, "model_params")
  expect_length(zlmp, 0)
  expect_identical(names(zlmp), character(0))
})

test_that("default model params generates correctly", {
  dparams <- default_model_params()
  expect_is(dparams, "model_params")

  init_names <- c("Pinit", "Ninit", "Zinit", "Dsinit", "kw", "kc", "tstep", "nyears",
                  "flag_stn", "flag_LI", "flag_atten", "flag_irrad", "flag_PIcurve",
                  "flag_grazing", "flag_outtype", "flag_outfreq", "flag_integ")
  parms_names <- c("VPmax", "alpha", "kN", "mP", "mp2", "Imax", "kz", "phiP",
                   "phiD", "betaz", "kNz", "mz", "mz2", "VD", "mD", "wmix", "CtoChl")

  expect_true(all(c(init_names, parms_names) %in% names(dparams)))
})

test_that("model params inherit correctly", {
  mp1 <- model_params(Pinit = 1, Ninit = 2, Zinit = 3)
  mp2 <- model_params(mp1, Pinit = "one", Ninit = "two")
  expect_equal(mp2$Pinit, "one")
  expect_equal(mp2$Ninit, "two")
  expect_equal(mp2$Zinit, 3)

  expect_identical(model_params(mp1), mp1)
})

test_that("model params can't inherit from anything other than a model_params", {
  expect_error(model_params(list()), "model_params must inherit from another model_params")
})
