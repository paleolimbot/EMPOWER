
# ------------------------------------------------------------------------------------- #
# L_I: Anderson (1993) function for calculating daily photosynthesis                    #
# Uses in-built piecewise light extinction and assumes sinusoidal irradiance during day #
# Water column divided into three layers: 0-5m, 5-23m and >23m                          #
# ------------------------------------------------------------------------------------- #

# FNaphy calculates values of a (chlorophyll abosorption) for each layer; alpha[layer] is then calculated as a[layer]*alphamax
# This function is only used when Anderson (1993) is selected for calculating daily depth-integrated photosynthesis

#' Anderson (1993) function for calculating daily photosynthesis
#'
#' Uses in-built piecewise light extinction and assumes sinusoidal irradiance during day
#' Water column divided into three layers: 0-5m, 5-23m and >23m
#'
#' @param ss Not sure
#' @param ztop not sure
#' @param zbottom not sure
#' @param aphylast not sure
#'
#' @return "acalc", whatever that is
#' @keywords internal
#'
FNaphy <- function(ss,ztop,zbottom,aphylast) {

  g <- c(0.048014,0.00023779,-0.023074,0.0031095,-0.0090545,0.0027974,0.00085217,-3.9804E-06,0.0012398,-0.00061991)   # coeffs. for calculating a#

  x <- zbottom+1.0
  xlg <- log(x)
  termf1a <- x*xlg-x
  termf2a <- x*xlg^2-2*x*xlg+2*x
  termf3a <- x*xlg^3-3*x*xlg^2+6*x*xlg-6*x
  x <- ztop+1.0
  xlg <- log(x)
  termf1b <- x*xlg-x
  termf2b <- x*xlg^2-2*x*xlg+2*x
  termf3b <- x*xlg^3-3*x*xlg^2+6*x*xlg-6*x

  terma <- g[1]+g[2]*ss+g[5]*ss^2+g[7]*ss^3
  termb <- g[3]+g[4]*ss+g[9]*ss^2
  termc <- g[6]+g[10]*ss
  acalc <- (zbottom+1.0)*terma+termf1a*termb+
    termf2a*termc+termf3a*g[8]-((ztop+1.0)*terma+
                                  termf1b*termb+termf2b*termc+termf3b*g[8])

  return(acalc)
}
