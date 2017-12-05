context("EMPOWER")

test_that("EMPOWER() default run generates original default output", {
  out_statevars <- read.csv(system.file("default_output/out_statevars.txt", package = "EMPOWER"),
                            header = FALSE)

  # EMPOWER() needs the working directory set specifically...
  # create a temporary directory in which to test
  tmp_dir <- "tmpdir"
  dir.create(tmp_dir)

  out2_rmse <- try(withr::with_dir(tmp_dir, {
    # copy default input files to temp directory
    input_files <- list.files(system.file("default_input", package = "EMPOWER"), full.names = TRUE)
    file.copy(input_files, basename(input_files))

    # run the model
    FNmodel_core()

    # read output file
    out2_statevars <- read.csv("out_statevars.txt", header = FALSE)

    # calculate RMSE, column wise
    sqrt(colMeans((out2_statevars - out_statevars)^2))
  }))

  # remove temporary directory
  unlink(tmp_dir, recursive = TRUE)

  expect_false(inherits(out2_rmse, "try-error"))

  # make sure all RMSE are less than 0.02
  expect_true(all(out2_rmse < 0.02))
})
