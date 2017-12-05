context("EMPOWER")

test_that("EMPOWER() with default params returns identical values as original model", {

  out_statevars <- read.csv(system.file("default_output/out_statevars.txt", package = "EMPOWER"),
                            header = FALSE)

  # run model with original default parameters
  empower_result <- EMPOWER(params = model_params(VPmax = 2.5, alpha = 0.15, kN = 0.85,
                                                  mP = 0.015, mp2 = 0.025, Imax = 1.0,
                                                  kz = 0.6, phiP = 0.67, phiD = 0.33,
                                                  betaz = 0.69, kNz = 0.75, mz = 0.02,
                                                  mz2 = 0.34, VD = 6.43, mD = 0.06,
                                                  wmix = 0.13, CtoChl = 75.0,
                                                  # these are 'init' values:
                                                  Pinit = 0.1, Ninit = 10, Zinit = 1.0,
                                                  Dsinit = 0.1, kw = 0.04, kc = 0.03,
                                                  tstep = 0.1, nyears = 5, flag_stn = 2,
                                                  flag_LI = 1, flag_atten = 2, flag_irrad = 2,
                                                  flag_PIcurve = 1, flag_grazing = 2,
                                                  flag_outtype = 1, flag_outfreq = 1,
                                                  flag_integ = 1))

  # extract statevars
  out2_statevars <- empower_result$statevars

  # calculate RMSE, column wise
  out2_rmse <- sqrt(colMeans((out2_statevars - out_statevars)^2))

  # make sure all RMSE are less than 0.02
  expect_true(all(out2_rmse < 0.02))
})
