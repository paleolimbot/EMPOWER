context("EMPOWER")

test_that("EMPOWER() default run generates original default output", {
  out_statevars <- read.csv(system.file("default_output/out_statevars.txt", package = "EMPOWER"),
                            header = FALSE)
  out_aux <- read.csv(system.file("default_output/out_aux.txt", package = "EMPOWER"),
                      header = FALSE)
  out_fluxes <- read.csv(system.file("default_output/out_fluxes.txt", package = "EMPOWER"),
                         header = FALSE)

  # currently EMPOWER() needs the working directory set...
  tmp_dir <- "tmpdir"
  dir.create(tmp_dir)
  old_wd <- getwd()
  setwd(tmp_dir)

  # copy default input files to temp directory
  input_files <- list.files(system.file("default_input", package = "EMPOWER"), full.names = TRUE)
  file.copy(input_files, basename(input_files))

  # run the model
  EMPOWER()

  # read output files
  out2_statevars <- read.csv("out_statevars.txt", header = FALSE)
  out2_aux <- read.csv("out_aux.txt", header = FALSE)
  out2_fluxes <- read.csv("out_fluxes.txt", header = FALSE)

  # expect equal to other outs
  expect_equal(out2_statevars, out_statevars)
  expect_equal(out2_aux, out_aux)
  expect_equal(out2_fluxes, out_fluxes)

  setwd(old_wd)
  unlink(tmp_dir, recursive = TRUE)
})
