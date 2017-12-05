
# ---------------------------------------------------------------------------------- #
# Function get_flux: the specification of the ecosystem model is primarily handled   #
#  here. The terms in the differential equations are calculated and transferred to   #
#  matrix, flux[i,j], which is passed back to the core code for integration.         #
#  The user can define auxiliary variables, Y[i], that are written to output         #
#  files (state variables are automatically stored for output).                      #
# ---------------------------------------------------------------------------------- #

#' Handle specification of the ecosystem model
#'
#' The terms in the differential equations are calculated and transferred to
#' matrix, flux[i,j], which is passed back to the core code for integration.
#' The user can define auxiliary variables, Y[i], that are written to output
#' files (state variables are automatically stored for output).
#'
#' @param X The array that carries state variables flux[i,j], which is passed back to the
#'   core code for integration
#'
#' @return A matrix flux[i,j], which is passed back to the core code for integration.
#' @keywords internal
#'
#' @references
#' Anderson (1993)
#' Evans and Parslow (1985)
#' Eppley function of sea surface temperature
#'
FNget_flux <- function(X) {      # X is the array that carries state variables
  # (note that all functions begin with "FN...", although this is not obligatory in R

  P <- X[1];N <- X[2];Z <- X[3]; Ds <- X[4]         # unpack state variables: P, N, Z, D

  # Environmental forcing for current day of year

  noonparnow <- FNnoonparcalc(daynow,latradians,clouds,e0)  # noon irradiance (PAR), W m-2
  daylnow <- FNdaylcalc(daynow,latradians)                  # day length, h
  MLD1 <- MLD[daynow]                                       # mixed layer depth at start of day
  MLD2 <- MLD[daynow+1]                                     # mixed layer depth at end of day
  MLDnow <- MLD1+(istepnow-1)/nstepsday*(MLD2-MLD1)         # MLD at start of current time step
  N0 <- aN*MLDnow + bN                                      # DIN immediately below mixed layer

  VpT <- Vp0*1.066^SST[daynow]             # maximum photosynthetic (growth) rate (Eppley function of sea surface temperature)

  ### Light attenuation in water column
  # Major choice here regarding calculation of light attenuation in the water column:
  # (1) Mixed layer as single layer, with light attenuation calculated according to k_tot = k_w + k_c*P
  # (2) Mixed layer split into three layers according to Anderson (1993): 0-5, 5-23, >23m, with a
  #       separate extinction coefficient for each layer as a function of chlorophyll

  chl <- P*6.625*12.0/CtoChl     # chlorophyll, mg m-3 (Redfield ratio of 6.625 mol C mol N-1 assumed for C:N of phytoplankton)
  ss <- sqrt(chl)                # square root of chlorophyll

  # Four unit array for the following variables handles three layers; first unit in array is for surface properties,
  #  elements 2,3,4 handle layers 1,2,3 respectively (if only one layer needed (option 1 or light attenuation), then
  #  elements 3,4 are redundant

  kPAR <- array(dim=c(4))        # light extinction coeff., m-1
  zbase <- array(dim=c(4))       # depth of base of layer; for element 1 in array (surface) this is zero
  Ibase <- array(dim=c(4))       # irradiance leaving the base of layer; for element 1, this is surface irradiance, W m-2
  zdep <- array(dim=c(4))        # depth of the layer, m
  aphybase <- array(dim=c(4))    # a (absorption by chlorophyll) at base of layer; only required when using Anderson (1993) calculation of photosynthesis

  # flag_atten is choice of light attenuation scheme: (1) single layer k_w+k_c*P, (2) three layer, extinction coefficients as fn chl (Anderson, 1993)

  if (flag_atten==1) {    # option (1)
    jnlay <- 1            # mixed layer as a single layer
    zbase[1] <- 0.
    zbase[2] <- MLDnow
    zdep[2] <- MLDnow
    kPAR[2] <- kw + kc*P
    Ibase[1] <- noonparnow                      # irradiance entering water column
    Ibase[2] <- Ibase[1]*exp(-kPAR[2]*MLDnow)   # irradiance at base of mixed layer

  } else if (flag_atten==2) {     # option (2): MLD separated into depth ranges 0-5m, 5-23m and >23m  Anderson (1993)

    Ibase[1] <- noonparnow   # irradiance entering water column

    kPAR[2] <- 0.13096 + 0.030969*ss + 0.042644*ss^2 - 0.013738*ss^3 + 0.0024617*ss^4 - 0.00018059*ss^5     # extinction coefficients (Anderson, 1993)
    kPAR[3] <- 0.041025 + 0.036211*ss + 0.062297*ss^2 - 0.030098*ss^3 + 0.0062597*ss^4 - 0.00051944*ss^5
    kPAR[4] <- 0.021517 + 0.050150*ss + 0.058900*ss^2 - 0.040539*ss^3 + 0.0087586*ss^4 - 0.00049476*ss^5

    zbase[1] <- 0.
    Ibase[1] <- noonparnow      # irradiance entering water column

    # Three layers only if MLD > 23.0m, otherwise one or two layers:

    if (MLDnow<=5.0) {
      jnlay <- 1
      zbase[2] <- MLDnow
      zdep[2] <- MLDnow
      Ibase[2] <- Ibase[1]*exp(-kPAR[2]*zdep[2])    # irradiance leaving layer 1
    } else if (MLDnow>5 && MLDnow<=23.0) {
      jnlay <- 2
      zbase[2] <- 5.0
      zdep[2] <- 5.0
      Ibase[2] <- Ibase[1]*exp(-kPAR[2]*5.0)        # irradiance leaving layer 1
      zbase[3] <- MLDnow
      zdep[3] <- MLDnow-5.0
      Ibase[3] <- Ibase[2]*exp(-kPAR[3]*zdep[3])    # irradiance leaving layer 2
    } else if (MLDnow>23.0) {
      jnlay <- 3
      zbase[2] <- 5.0
      zdep[2] <- 5.0
      Ibase[2] <- Ibase[1]*exp(-kPAR[2]*5.0)        # irradiance leaving layer 1
      zbase[3] <- 23.0
      zdep[3] <- 23.0-5.0
      Ibase[3] <- Ibase[2]*exp(-kPAR[3]*zdep[3])    # irradiance leaving layer 2
      zbase[4] <- MLDnow
      zdep[4] <- MLDnow-23.0
      Ibase[4] <- Ibase[3]*exp(-kPAR[4]*zdep[4])    # irradiance leaving layer 3
    }

  } # end if: flag_atten

  ### Calculate L_I (light limitation of growth, 0 <= L_I <= 1)

  L_Isum <- 0.     # L_I is calculated as a weighted sum over the total mixed layer

  for (ilay in seq(2,jnlay+1)) {         # loop over layers; element 2 in array corresponds to first layer

    # Call function for calculating photosynthesis

    if (flag_LI==1) {                      # numeric integration for light over time (through day) with analytic depth integrals
      L_I <- FNLIcalcNum(zdep[ilay],Ibase[ilay-1],Ibase[ilay],kPAR[ilay],alpha,VpT,daylnow,flag_irrad,flag_PIcurve)

    } else if (flag_LI==2) {               # Evans and Parslow (1985): triangular light, Smith fn for P-I curve
      L_I <- FNLIcalcEP85(zdep[ilay],Ibase[ilay-1],Ibase[ilay],kPAR[ilay],alpha,VpT,daylnow)

    }  else if (flag_LI==3) {              # Anderson (1993): sinusoidal light, exponential fn for P-I curve, alpha spectrally dependent
      aphybase[1] <- 0.36796 + 0.17537*ss - 0.065276*ss^2 + 0.013528*ss^3 - 0.0011108*ss^4           # a (chl absorption) at ocean surface as function chl
      ahash <- FNaphy(ss,zbase[ilay-1],zbase[ilay],aphybase[ilay-1])                                 # change in a with depth
      aphybase[ilay] <- aphybase[ilay-1]+ahash                                                       # a at base of layer
      aphyav <- aphybase[ilay-1]+ahash*0.5                                                           # average a in layer (from which alpha is calculated: alpha = a*alphamax
      L_I <- FNLIcalcA93(zdep[ilay],Ibase[ilay-1],Ibase[ilay],kPAR[ilay],alpha,VpT,daylnow,aphyav)

    }

    L_Isum <- L_Isum + L_I*zdep[ilay]      # multiply by layer depth in order to set up weighted average for total mixed layer

  } # end loop: layers

  L_I <- L_Isum/MLDnow                   # weighted average for mixed layer

  ### Calculate L_N (nutrient limitation of growth, 0 <= L_I <= 1)

  L_N <- N/(kN+N)                  # nutrient limitation of growth rate (0 <= L_N <= 1)

  ### Grazing

  if (flag_grazing==1) {           # phihat_i = phi_i (grazing preferences)
    phiPhat <- phiP
    phiDhat <- phiD
  } else {
    phiPhat <- phiP*P              # phihat_i = phi_i*P_i
    phiDhat <- phiD*Ds
  }
  intakespP <- Imax*phiPhat*P/(kz^2+phiPhat*P+phiDhat*Ds)    # specific intake: phytoplankton, d-1
  intakespD <- Imax*phiDhat*Ds/(kz^2+phiPhat*P+phiDhat*Ds)   # specific intake: detritus, d-1

  # terms in the differential equations  (unit: mmol N m-3 d-1)

  Pgrow <- VpT*L_I*L_N*24/CtoChl*P           # P growth
  Pgraz <- intakespP*Z                       # P grazed
  Dgraz <- intakespD*Z                       # D grazed
  Pmort <- mP*P                              # P mortality: linear
  Pmort2 <- mP2*P*P                          # P mortality: quadratic
  Pmix  <- wmix*P/MLDnow                     # P loss: mixing
  Zmix  <- wmix*Z/MLDnow                     # Z loss: mixing
  Nmix  <- wmix*(N0-N)/MLDnow                # N net input: mixing
  Zgrow <- betaz*kNz*(Pgraz+Dgraz)           # Z growth
  Zexc  <- betaz*(1.0-kNz)*(Pgraz+Dgraz)     # Z excretion
  Zpel  <- (1.0-betaz)*(Pgraz+Dgraz)         # Z faecal pettet production
  Zmort <- mz*Z                              # Z mortality: linear
  Zmort2 <- mz2*Z^2                          # Z mortality: quadratic
  Dmix  <- wmix*Ds/MLDnow                    # D loss: mixing
  Dsink <- VD*Ds/MLDnow                      # D loss: sinking
  Dremin <- mD*Ds                            # D loss: remineralisation

  #### Entrainment/dilution (when depth of mixed layer increases)

  if (MLD2-MLD1>0.0) {
    Nadd <- (MLD2-MLD1)*N0/MLDnow       # entrainment: N
    Ndilute <- N*(MLD2-MLD1)/MLDnow     # dilution: N already present in ML
    Pdilute <- P*(MLD2-MLD1)/MLDnow     # dilution: P
    Zdilute <- Z*(MLD2-MLD1)/MLDnow     # dilution: Z
    Ddilute <- Ds*(MLD2-MLD1)/MLDnow    # dilution: Ds
  } else {
    Nadd <- 0.0                         # detrainment, but concentration in ML unchanged
    Ndilute <- 0.0
    Pdilute <- 0.0
    Zdilute <- 0.0
    Ddilute <- 0.0
  }

  ### Transfer terms as calculated above are trasferred to matrix flux[i,j] for processing by the core code
  # (number of terms for any one state variable should not exceed nfluxmax as specified in core code)

  # Phytoplankton
  flux[1,1] <- Pgrow                 # P growth
  flux[2,1] <- -Pgraz                # P grazed
  flux[3,1] <- -Pmort                # P non-grazing mortality: linear
  flux[4,1] <- -Pmort2               # P non-grazing mortality: quadratic
  flux[5,1] <- -Pmix                 # P loss: mixing
  flux[6,1] <- -Pdilute              # P loss: dilution


  # Nitrate
  flux[1,2] <- Nmix                  # N net input: mixing
  flux[2,2] <- Zexc                  # Z excretion
  flux[3,2] <- Dremin                # D remineralisation
  flux[4,2] <- -Pgrow                # uptake by P
  flux[5,2] <- Nadd                  # N input: entrainment
  flux[6,2] <- -Ndilute              # N loss: dilution

  # Zooplankton
  flux[1,3] <- Zgrow                 # Z growth
  flux[2,3] <- -Zmort                # Z mortality: linear
  flux[3,3] <- -Zmort2               # Z mortality: quadratic
  flux[4,3] <- -Zdilute              # Z loss: dilution
  flux[5,3] <- -Zmix                 # Z loss: mixing

  # Detritus
  flux[1,4] <- Pmort                 # D from P mortality (linear)
  flux[2,4] <- Pmort2                # D from P mortality (quadratic)
  flux[3,4] <- Zmort                 # D from Z mortality (linear)
  flux[4,4] <- Zpel                  # D from Z faecal pellets
  flux[5,4] <- -Dgraz                # D grazed
  flux[6,4] <- -Dmix                 # D loss: mixing
  flux[7,4] <- -Dsink                # D loss: sinking
  flux[8,4] <- -Dremin               # D loss: remineralisation
  flux[9,4] <- -Ddilute              # D loss: dilution

  # User-definied auxiliary variables for writing to output files: number stored = nDvar (core code)
  # If you want to increase or decrease the number of units in this array, this is controlled by nDvar in core code

  Y[1]  <<- L_I
  Y[2]  <<- L_N
  Y[3]  <<- Pgrow/P
  Y[4]  <<- MLDnow
  Y[5]  <<- noonparnow
  Y[6]  <<- chl
  Y[7]  <<- SST[daynow]
  Y[8]  <<- Pgrow*MLDnow
  Y[9]  <<- Pgrow
  Y[10] <<- Pgraz
  Y[11] <<- Pmort+Pmort2
  Y[12] <<- 0
  Y[13] <<- 0
  Y[14] <<- 0
  Y[15] <<- 0
  Y[16] <<- 0
  Y[17] <<- 0
  Y[18] <<- 0
  Y[19] <<- 0
  Y[20] <<- 0

  return(flux)      # matrix flux is passed back to the core code
}
