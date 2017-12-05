
#' Create model parameter specifications
#'
#' Creates a default parameter specification as distributed in the original file
#' "NPZD_parms.txt" and "NPZD_extras.txt". Values can be overridden from default
#' using default_model_params(), or overridden from an arbitrary parameter
#' specification using model_params().
#'
#' @param VPmax VPmax(0) (2.5): max. rate photosynthesis, g C (g chl)-1 h-1
#' @param alpha alpha (0.15): initial slope of P-I curve, g C (g chl)-1 h-1 (W
#'   m-2)-1
#' @param kN kN (0.85): half saturation constant for N uptake, mmol m-3
#' @param mP mP (0.02): phytoplankton mortality rate, linear d-1
#' @param mp2 mp2 (0.025): phytoplankton mortality rate, quadratic, (mmol N
#'   m-3)-1 d-1
#' @param Imax Imax (1.0): zooplankton max. ingestion, d-1
#' @param kz kz (0.86): zooplankton half sat. for intake, mmol N m-3
#' @param phiP phiP (0.67): zooplankton preference: P
#' @param phiD phiD (0.33): zooplankton preference: D
#' @param betaz betaz (0.69): zooplankton absorption efficiency: N
#' @param kNz kNz (0.75): zooplankton net production efficiency: N
#' @param mz mz (0.02): zooplankton mortality rate, linear, d-1
#' @param mz2 mz2 (0.34): zooplankton moratlity, quadratic, (mmol N m-3)-1 d-1
#' @param VD VD (6.43): detritus sinking rate, m d-1
#' @param mD mD (0.06): deteritus remineralisation rate, d-1
#' @param wmix wmix (0.13): cross-thermocline mixing rate, m d-1
#' @param CtoChl CtoChl (75.0) carbon to chlorophyll ratio, g g-1
#'
#' @param Pinit initial P, mmol N m-3
#' @param Ninit initial N, mmol m-3
#' @param Zinit initial Z, mmol N m-3
#' @param Dsinit inital Ds, mmol N m-3
#' @param kw light attenuation by water, m-1
#' @param kc light attenuation by phyto., (mmol N m-3)-1
#' @param tstep time step, day
#' @param nyears run duration, years (1 year= 365 days)
#' @param flag_stn station; 1=India (60N 20W), 2=Biotrans (47N 20W), 3=Kerfix
#'   (50 40S 68 25E), 4=Papa (50N, 145W)
#' @param flag_LI 1=numeric; 2=Evans & Parslow (1985); 3=Anderson (1993)
#' @param flag_atten 1: kw+kcP (one layer for ML); 2: Anderson (1993) (spectral,
#'   layers 0-5,5-23,>23 m
#' @param flag_irrad 1=triangular over day; 2=sinusoidal over day
#' @param flag_PIcurve 1=Smith fn; 2=exponential fn
#' @param flag_grazing (1) phihat_i = phi_i, (2) phihat_i = phi_i*P_i
#' @param flag_outtype output to files: 0=none, 1=last year only, 2=whole
#'   simulation
#' @param flag_outfreq frequency for writing output: 1=every day, 0=every time
#'   step
#' @param flag_integ choice of integration method: 0=Euler, 1=Runge Kutta 4
#'
#' @return A model_params object, which is a list of parameter values.
#' @export
#'
#' @examples
#' default_model_params()
#'
default_model_params <- function() {

  model_params(.inherit_from = NULL,
               VPmax = 2.5, alpha = 0.15, kN = 0.85,
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
               flag_integ = 1)
}

#' @rdname default_model_params
#' @export
model_params <- function(.inherit_from = NULL,
                         # these are 'parms' values
                         VPmax = NULL, alpha = NULL, kN = NULL,
                         mP = NULL, mp2 = NULL, Imax = NULL,
                         kz = NULL, phiP = NULL, phiD = NULL,
                         betaz = NULL, kNz = NULL, mz = NULL,
                         mz2 = NULL, VD = NULL, mD = NULL,
                         wmix = NULL, CtoChl = NULL,
                         # these are 'init' values:
                         Pinit = NULL, Ninit = NULL, Zinit = NULL,
                         Dsinit = NULL, kw = NULL, kc = NULL,
                         tstep = NULL, nyears = NULL, flag_stn = NULL,
                         flag_LI = NULL, flag_atten = NULL, flag_irrad = NULL,
                         flag_PIcurve = NULL, flag_grazing = NULL,
                         flag_outtype = NULL, flag_outfreq = NULL,
                         flag_integ = NULL) {

  if(is.null(.inherit_from)) {
    .inherit_from <- list()
  } else if(!inherits(.inherit_from, "model_params")) {
    stop("model_params must inherit from another model_params")
  }

  # create list() of new params
  new_params <-  list(VPmax = VPmax, alpha = alpha, kN = kN,
                      mP = mP, mp2 = mp2, Imax = Imax,
                      kz = kz, phiP = phiP, phiD = phiD,
                      betaz = betaz, kNz = kNz, mz = mz,
                      mz2 = mz2, VD = VD, mD = mD,
                      wmix = wmix, CtoChl = CtoChl,
                      # these are 'init' values:
                      Pinit = Pinit, Ninit = Ninit, Zinit = Zinit,
                      Dsinit = Dsinit, kw = kw, kc = kc,
                      tstep = tstep, nyears = nyears, flag_stn = flag_stn,
                      flag_LI = flag_LI, flag_atten = flag_atten, flag_irrad = flag_irrad,
                      flag_PIcurve = flag_PIcurve, flag_grazing = flag_grazing,
                      flag_outtype = flag_outtype, flag_outfreq = flag_outfreq,
                      flag_integ = flag_integ)

  # remove NULLs
  new_params <- new_params[!vapply(new_params, is.null, logical(1))]

  # check for names
  if(length(new_params) != 0 && is.null(names(new_params)) || any(names(new_params) == "")) {
    stop("model_params must have all named arguments")
  }

  # combine with previous list, subset by unique name to perform override
  new_params <- c(new_params, .inherit_from)
  new_params <- new_params[unique(names(new_params))]

  # return list with a custom class
  new_params <- structure(new_params, class = "model_params", .Names = as.character(names(new_params)))

  # return model_params
  return(new_params)
}

#' Print a model_params object
#'
#' @param x An object
#' @param ... ignored
#'
#' @return The object, invisibly
#' @export
#'
print.model_params <- function(x, ...) {
  cat(sprintf("<%s>\n", class(x)[1]))
  if(length(x) == 0) {
    cat("<empty>")
  } else {
    print(as.data.frame(unclass(x)))
  }
  invisible(x)
}
