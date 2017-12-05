
# ----------------------------------------------------------------------- #
# Calculation of day length as function of day of year and latitude       #
# ----------------------------------------------------------------------- #

#' Calculation of day length as function of day of year and latitude
#'
#' @param jday The julian day of the year
#' @param latradians Latitude, in radians
#'
#' @return Day length (in hours?)
#' @keywords internal
#'
FNdaylcalc <- function(jday,latradians) {
  declin <- 23.45*sin(2*pi*(284+jday)*0.00274)*pi/180      # solar declination angle
  daylnow <- 2*acos(-1*tan(latradians)*tan(declin))*12/pi
  return(daylnow)
}
