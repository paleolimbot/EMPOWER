
# -------------------------------------------------- #
# Calculation of noon irradiance                     #
# -------------------------------------------------- #

# Noon PAR, W m-2

#' Calculation of noon irradiance
#'
#' @param jday Julian day of the year
#' @param latradians Latitude (in radians)
#' @param clouds Clouds parameter
#' @param e0 no idea
#'
#' @return Noon PAR, W m-2
#' @keywords internal
#'
FNnoonparcalc <- function(jday,latradians,clouds,e0) {
  albedo <- 0.04                # albedo
  solarconst <- 1368.0          # solar constant, w m-2
  parrac <- 0.43                # PAR fraction
  declin <- 23.45*sin(2*pi*(284+jday)*0.00274)*pi/180                # solar declination angle
  coszen <- sin(latradians)*sin(declin)+cos(latradians)*cos(declin)  # cosine of zenith angle
  zen <- acos(coszen)*180/pi                                         # zenith angle, degrees
  Rvector <- 1/sqrt(1+0.033*cos(2*pi*jday*0.00274))                  # Earth's radius vector
  Iclear <- solarconst*coszen^2/(Rvector^2)/(1.2*coszen+e0*(1.0+coszen)*0.001+0.0455)   # irradiance at ocean surface, clear sky
  cfac <- (1-0.62*clouds*0.125+0.0019*(90-zen))                      # cloud factor (atmospheric transmission)
  Inoon <- Iclear*cfac*(1-albedo)                                    # noon irradiance: total solar
  noonparnow <- parrac*Inoon

  return(noonparnow)
}
