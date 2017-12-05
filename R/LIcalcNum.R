
# ------------------------------------------------------------------------------------------------------ #
# L_I: Calculation of daily photosynthesis: numeric intergration (time) and analytic integration (depth) #
# Diel cycle of irradiance: choice of (1) triangular vs (2) sinusoidal (flag_irrad)                      #
# P-I curve: choice of (1) Smith fn vs (2) exponential fn  (flag_PIcurve)                                #
# numeric integration (time) uses nstepst+1 steps and trapezoidal rule for first and last steps          #
#   (the user may manually change nstepst in the code if so desired)                                     #
# -------------------------------------------------------------------------------------------------------#

#' Calculation of daily photosynthesis
#'
#' numeric intergration (time) and analytic integration (depth)
#' Diel cycle of irradiance: choice of (1) triangular vs (2) sinusoidal (flag_irrad)
#' P-I curve: choice of (1) Smith fn vs (2) exponential fn  (flag_PIcurve)
#' numeric integration (time) uses nstepst+1 steps and trapezoidal rule for first and last steps
#' (the user may manually change nstepst in the code if so desired)
#'
#' @param zdepth no idea
#' @param Iin no idea
#' @param Iout no idea
#' @param kPARlay no idea
#' @param alpha no idea
#' @param Vp no idea
#' @param daylnow no idea
#' @param choiceirrad no idea
#' @param choicePI no idea
#'
#' @return Lim_I (probably a scalar?)
#' @keywords internal
#'
FNLIcalcNum <- function(zdepth,Iin,Iout,kPARlay,alpha,Vp,daylnow,choiceirrad,choicePI) {

  sumps <- 0.0                     # for calculating sum of photosynthesis
  range_rad = pi/2.0
  nstepst <- 10                    # the day from sunrise to midday will be divided into nstepst steps (can be altered if so desired)

  for (itme in seq(1,nstepst)) {   # time loop     (no need to calculate itme = 0 because irradiance zero at sunrise)

    if(itme==nstepst) {            # last step (note use of trapezoidal rule for integration)
      rtrapt <- 0.5 } else {
        rtrapt <- 1.0                # other steps
      }

    if (choiceirrad==1) {            # coice of either (1) triangular or (2) sinusoidal irradiance over day
      Icorr <- itme/nstepst                    # irradiance multiplier for time of day: triangular
    } else if (choiceirrad==2) {
      Icorr <- sin(itme/nstepst*range_rad)     # irradiance multiplier for time of day: sinusoidal
    }

    if (choicePI==1) {               # Smith function (the most straightforward choice in terms of calculating analytic depth integral)

      sumpsz <- 0.
      x0 <- alpha*Icorr*Iin
      xH <- alpha*Icorr*Iin*exp(-kPARlay*zdepth)
      sumpsz <- Vp/kPARlay/zdepth*(log(x0+(Vp^2+x0^2)^0.5)-log(xH+(Vp^2+xH^2)^0.5))

    } else {                        # exponential function (Platt et al., 1990)

      # Default setting here is analytic depth integral. Based on a factorial sequence, we found that this does not always perform well and so, if desired,
      # the user can switch to an analytic depth integral by setting (manually in the code) choiceanalytic <- 2

      choiceanalytic <- 1         # 1=analytic solution; 2=numeric solution; set manually in the code

      if (choiceanalytic==1) {    # analytic solution
        sumpsz <- 0.
        istar1 <- alpha*Icorr*Iin/Vp
        istar2 <- alpha*Icorr*Iout/Vp

        # Analytic solution is an infinite series
        for (ifact in seq(1,16)) {  # 16 units in factorial array (by all means increase this number by altering the code; run time will increase accordingly)
          sumpsz <- sumpsz + (-1)^(ifact+1)/(ifact*factorial(ifact))*(istar1^ifact-istar2^ifact)
        }

        sumpsz <- sumpsz*Vp/kPARlay/zdepth

      } else {   # numeric solution

        nstepsz <- 100            #no. of steps: depth
        sumpsz <- 0.0

        for (izz in seq(0,nstepsz)) {     # depth loop
          if(izz==0 | izz==nstepsz) {
            rtrapz <- 0.5 } else {        # trapezoidal rule
              rtrapz <- 1.0
            }

          z <- izz/nstepsz*zdepth
          Iz <- Iin*Icorr*exp(-kPARlay*z)

          psnow <- Vp*(1.0-exp(-alpha*Iz/Vp))

          sumpsz <- sumpsz + psnow*rtrapz

        }  # end loop: depth

        sumpsz <- sumpsz/nstepsz

      } # choice analytic/numeric

    } # end if: choice PI curve

    sumps <- sumps + sumpsz*rtrapt       # sum photosynthesis over daylight hours

  }  # end loop: time

  Lim_I <- sumps/nstepst*daylnow/24.0  # take account of zero PS during night
  Lim_I <- Lim_I/Vp                    # dimensionless
  return(Lim_I)
}
