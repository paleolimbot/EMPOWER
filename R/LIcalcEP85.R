
# ------------------------------------------------------------------------------ #
# L_I: Evans & Parslow (1985) function for calculating daily photosynthesis      #
# Light extinction using Beer's Law and assumed triangular irradiance during day #
# ------------------------------------------------------------------------------ #

#' Evans & Parslow (1985) function for calculating daily photosynthesis
#'
#' Light extinction using Beer's Law and assumed triangular irradiance during day
#'
#' @param zdepth not sure
#' @param Iin not sure
#' @param Iout not sure
#' @param kPARlay not sure
#' @param alpha not sure
#' @param Vp not sure
#' @param daylnow not sure
#'
#' @return Lim_I (probably a scalar)
#' @keywords internal
#'
#' @references
#' Evans & Parslow (1985)
#'
FNLIcalcEP85 <- function(zdepth,Iin,Iout,kPARlay,alpha,Vp,daylnow) {
  tau1 <- daylnow/2.0/24.0
  betaEP <- Vp*tau1/alpha/Iin
  betaEP2 <- Vp*tau1/alpha/Iout
  Lim_I <- 2.0*Vp/zdepth/kPARlay*(sqrt(betaEP2^2+tau1^2) -
    tau1*log((tau1+sqrt(betaEP2^2+tau1^2))/betaEP2) -
    (sqrt(betaEP^2+tau1^2)-tau1*log((tau1+sqrt(betaEP^2+tau1^2))/betaEP)) -
      betaEP2+betaEP)

  Lim_I <- Lim_I/Vp  # Divide by Vp to make dimensionless
  return(Lim_I)
}

