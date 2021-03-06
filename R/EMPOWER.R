
#' Run the EMPOWER model
#'
#' @param params A parameter set generated by \link{model_params}
#' @param directory An optional working directory in which to run the model
#'
#' @return An object of class EMPOWER_result
#' @export
#'
EMPOWER <- function(params = default_model_params(), directory = NULL) {

  # check params object
  null_args <- names(params)[vapply(params, is.null, logical(1))]
  if(length(null_args) > 0) stop("Correct the following NULL arguments: ", paste0(null_args, collapse = ", "))

  # if directory is NULL, create a temporary one and delete it on exit
  if(is.null(directory)) {
    directory <- tempfile()
    on.exit(unlink(directory, recursive = TRUE))
  }

  # try to create the output directory if it doesn't exist
  if(!dir.exists(directory)) dir.create(directory)
  if(!dir.exists(directory)) stop("Could not create directory: ", directory)

  # get list of parameter names for each input file
  init_names <- c("Pinit", "Ninit", "Zinit", "Dsinit", "kw", "kc", "tstep", "nyears",
                  "flag_stn", "flag_LI", "flag_atten", "flag_irrad", "flag_PIcurve",
                  "flag_grazing", "flag_outtype", "flag_outfreq", "flag_integ")
  parms_names <- c("VPmax", "alpha", "kN", "mP", "mp2", "Imax", "kz", "phiP",
                   "phiD", "betaz", "kNz", "mz", "mz2", "VD", "mD", "wmix", "CtoChl")

  # create data frames that will be turned into NPZD_parms.txt and NPZD_extra.txt
  init_df <- data.frame(value = unlist(params[init_names]),
                        parameter = init_names, stringsAsFactors = FALSE)
  parms_df <- data.frame(value = unlist(params[parms_names]),
                         parameter = parms_names, stringsAsFactors = FALSE)

  # switch the working directory for call to FNmodel_core()
  withr::with_dir(directory, {
    # write extra and parms files
    write.csv(init_df, "NPZD_extra.txt", row.names = FALSE)
    write.csv(parms_df, "NPZD_parms.txt", row.names = FALSE)

    # copy stations_forcing.txt to directory
    file.copy(system.file("default_input/stations_forcing.txt", package = "EMPOWER"), "stations_forcing.txt")

    # run FNmodel_core()
    FNmodel_core()

    # extract output files out_aux.txt, out_fluxes.txt, out_statevars.txt
    # (all headerless CSV files)
    aux <- tibble::as_tibble(read.csv("out_aux.txt", header = FALSE))
    fluxes <- tibble::as_tibble(read.csv("out_fluxes.txt", header = FALSE))
    statevars <- tibble::as_tibble(read.csv("out_statevars.txt", header = FALSE))

    # return all as a classed list
    structure(list(aux = aux, fluxes = fluxes, statevars = statevars, params = params),
              class = "EMPOWER_result")
  })
}
