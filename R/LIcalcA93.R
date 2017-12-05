
# Photosynthesis calculated using polynomial approximation (Anderson, 1993)
# This function is only used when Anderson (1993) is selected for calculating daily depth-integrated photosynthesis

#' Polynomial approximation of daily? phyotosynthesis
#'
#' Photosynthesis calculated using polynomial approximation (Anderson, 1993)
#' This function is only used when Anderson (1993) is selected for calculating daily depth-integrated photosynthesis
#'
#' @param zdepth no idea
#' @param Iin no idea
#' @param Iout no idea
#' @param kPARlay no idea
#' @param alpha no idea
#' @param Vp no idea
#' @param daylnow no idea
#' @param ahashnow no idea
#'
#' @return Lim_I, probably a scalar
#' @keywords internal
#'
#' @references
#' Anderson (1993)
#' Platt et al. (1990)
#'
FNLIcalcA93 <- function(zdepth,Iin,Iout,kPARlay,alpha,Vp,daylnow,ahashnow) {

  omeg <- c(1.9004,-2.8333E-01,2.8050E-02,-1.4729E-03,3.0841E-05)        # polynomial coefficients for calculating photosynthesis (Platt et al., 1990)

  # Calculate alphamax
  alphamax <- alpha*2.602     # alphamax is alpha at wavelength of maximum absorption cross section

  # Calculate daily photosynthesis in each layer, mg C m-2 d-1
  V0 <- daylnow*Vp/(pi*kPARlay)
  V1 <- alphamax*ahashnow*Iin/Vp
  V2 <- alphamax*ahashnow*Iout/Vp
  Qpsnow <- V0*(omeg[1]*(V1-V2)+omeg[2]*(V1^2-V2^2)+omeg[3]*(V1^3-V2^3)+omeg[4]*(V1^4-V2^4)+omeg[5]*(V1^5-V2^5))

  # convert to dimensionless units to get 0 <= J <= 1
  Lim_I <- Qpsnow/zdepth/24.0/Vp

  return(Lim_I)
}

