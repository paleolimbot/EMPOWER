###################################################################
# --------------------------------------------------------------- #
# EMPOWER-1.0: Tom Anderson 2015                                  #
# --------------------------------------------------------------- #
###################################################################

# ---------------------------------------------------------------------------------- #
# Note that R functions are listed before the core code.                             #
# Ensure that input files are within the same directory as the model code:           #
#   NPZDparms.txt, NPZDextra.txt, stations_forcing.txt                               #
# ---------------------------------------------------------------------------------- #

rm(list=ls()) #clear workspace

# ---------------------------------------------------------------------------------- #
# Function get_flux: the specification of the ecosystem model is primarily handled   #
#  here. The terms in the differential equations are calculated and transferred to   #
#  matrix, flux[i,j], which is passed back to the core code for integration.         #
#  The user can define auxiliary variables, Y[i], that are written to output         #
#  files (state variables are automatically stored for output).                      #
# ---------------------------------------------------------------------------------- #

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

# ------------------------------------------------------------------------------------------------------ #
# L_I: Calculation of daily photosynthesis: numeric intergration (time) and analytic integration (depth) #
# Diel cycle of irradiance: choice of (1) triangular vs (2) sinusoidal (flag_irrad)                      #
# P-I curve: choice of (1) Smith fn vs (2) exponential fn  (flag_PIcurve)                                #
# numeric integration (time) uses nstepst+1 steps and trapezoidal rule for first and last steps          #
#   (the user may manually change nstepst in the code if so desired)                                     #
# -------------------------------------------------------------------------------------------------------#

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


# ------------------------------------------------------------------------------ #
# L_I: Evans & Parslow (1985) function for calculating daily photosynthesis      #
# Light extinction using Beer's Law and assumed triangular irradiance during day #
# ------------------------------------------------------------------------------ #

FNLIcalcEP85 <- function(zdepth,Iin,Iout,kPARlay,alpha,Vp,daylnow) {
tau1 <- daylnow/2.0/24.0
betaEP <- Vp*tau1/alpha/Iin
betaEP2 <- Vp*tau1/alpha/Iout
Lim_I <- 2.0*Vp/zdepth/kPARlay*(sqrt(betaEP2^2+tau1^2)
+  -tau1*log((tau1+sqrt(betaEP2^2+tau1^2))/betaEP2)
+  -(sqrt(betaEP^2+tau1^2)-tau1*log((tau1+
+  sqrt(betaEP^2+tau1^2))/betaEP))-betaEP2+betaEP)

Lim_I <- Lim_I/Vp  # Divide by Vp to make dimensionless
return(Lim_I)
}

# ------------------------------------------------------------------------------------- #
# L_I: Anderson (1993) function for calculating daily photosynthesis                    #
# Uses in-built piecewise light extinction and assumes sinusoidal irradiance during day #
# Water column divided into three layers: 0-5m, 5-23m and >23m                          #
# ------------------------------------------------------------------------------------- #

# FNaphy calculates values of a (chlorophyll abosorption) for each layer; alpha[layer] is then calculated as a[layer]*alphamax
# This function is only used when Anderson (1993) is selected for calculating daily depth-integrated photosynthesis

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

# Photosynthesis calculated using polynomial approximation (Anderson, 1993)
# This function is only used when Anderson (1993) is selected for calculating daily depth-integrated photosynthesis

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

# ----------------------------------------------------------------------- #
# Calculation of day length as function of day of year and latitude       #
# ----------------------------------------------------------------------- #

FNdaylcalc <- function(jday,latradians) {
declin <- 23.45*sin(2*pi*(284+jday)*0.00274)*pi/180      # solar declination angle
daylnow <- 2*acos(-1*tan(latradians)*tan(declin))*12/pi
return(daylnow)
}

# -------------------------------------------------- #
# Calculation of noon irradiance                     #
# -------------------------------------------------- #

# Noon PAR, W m-2
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

# End of functions

#################################################
# --------------------------------------------- #
#         Core model code starts here           #
# --------------------------------------------- #
#################################################

# --------------------------------------------- #
# Read in parameter values from file            #
# --------------------------------------------- #

theta <- read.csv("NPZD_parms.txt", quote="'")

# Assign parameter names and make parameters universal in the code (<<-)
Vp0         <<- theta[1,1]     # max rate photosynthesis at zero degrees C, g C (g chl)-1 h-1
alpha       <<- theta[2,1]     # initial slope of P-I curve, g C (g chl)-1 h-1 (W m-2)-1
kN          <<- theta[3,1]     # half sat. constant: N, mmol m-3
mP          <<- theta[4,1]     # phyto. mortality rate, linear, d-1
mP2         <<- theta[5,1]     # phyto. mortality rate, quadratic, (mmol N m-3)-1 d-1
Imax        <<- theta[6,1]     # zooplankton max. ingestion, d-1
kz          <<- theta[7,1]     # zoo. half sat. for intake, mmol N m-3
phiP        <<- theta[8,1]     # zoo. preference: P
phiD        <<- theta[9,1]     # zoo. preference: D
betaz       <<- theta[10,1]    # zoo. absorption efficiency: N
kNz         <<- theta[11,1]    # zoo. net production efficiency: N
mz          <<- theta[12,1]    # zoo. mortality rate, linear, d-1
mz2         <<- theta[13,1]    # zoo. mortality, quadratic, (mmol N m-3)-1 d-1
VD          <<- theta[14,1]    # detritus sinking rate, m d-1
mD          <<- theta[15,1]    # detritus remineralisation rate, d-1
wmix        <<- theta[16,1]    # cross-thermocline mixing rate, m d-1
CtoChl      <<- theta[17,1]    # carbon to chlorophyll ratio, g g-1

# -------------------------------------------------- #
# Read in initial conditions and run characteristics #
# -------------------------------------------------- #
    
init <- read.csv("NPZD_extra.txt", quote="'")
Pinit       <<- init[1,1]     # initial P, mmol N m-3
Ninit       <<- init[2,1]     # initial N, mmol N m-3
Zinit       <<- init[3,1]     # inital  Z, mmol N m-3
Dsinit      <<- init[4,1]     # initial Ds, mmol N m-3
kw          <<- init[5,1]     # light attenuation coeff.: water,  m-1
kc          <<- init[6,1]     # light attenuation coeff.: phytoplankton, m2 (mmol N)-1
delta_t     <<- init[7,1]     # time step, d (choose at time step that divides exactly into 1 as an integer)
nyears      <<- init[8,1]     # run duration, years
flag_stn    <<- init[9,1]     # choice of station; 1=India (60N 20W), 2=Biotrans (47N 20W), 3=Kerfix (50 40’S 68 25’E), 4=Papa (50N, 145W)
flag_LI     <<- init[10,1]    # 1=numeric; 2=Evans & Parslow (1985); 3=Anderson (1993)
flag_atten  <<- init[11,1]    # 1: kw+kcP (one layer for ML); 2: Anderson (1993) (spectral, layers 0-5,5-23,>23 m
flag_irrad  <<- init[12,1]    # 1=triangular over day; 2=sinusoidal over day
flag_PIcurve <<- init[13,1]   # 1=Smith fn; 2=exponential fn
flag_grazing <<- init[14,1]   # (1) phihat_i = phi_i, (2) phihat_i = phi_i*P_i
flag_outtype <<- init[15,1]   # output choice: 0=none, 1=last year only, 2=whole simulation
flag_outfreq <<- init[16,1]   # output frequency: 1=every day, 0=every time step
flag_integ <<- init[17,1]     # choice of integration method: 0=Euler, 1=Runge Kutta 4
    
print("**************************************************")
print("                                                  ")
print("input: model parameters read in from file NPZD_parms.txt")
print("input: model run characteristics read in from file NPZDextra.txt")
print("                                                  ")

# --------------------------------------------- #
# Settings dependent on choice of station       #
# --------------------------------------------- #  

MLDi <- array(dim=c(13))        # array to hold monthly mixed layer depths (from data)
SSTi <- array(dim=c(13))        # array to hold monthly sea surface temperature (from data)

forcing <- read.csv("stations_forcing.txt", header=TRUE)      # forcing data from text file (comma separated)

if (flag_stn==1) {                  # Station India (60N 20W)
  for (i in seq(1,13)) {
  MLDi[i] <- forcing[i,1]      # i,1 is first column in the text file
  SSTi[i] <- forcing[i,2]      # i,2 is second column in the text file
  }
  latitude <<- 60.0            # latitude, degrees
  clouds <<- 6.0               # cloud fraction, oktas
  e0 <<- 12.0                  # atmospheric vapour pressure
  aN <<- 0.0074                # coeff. for N0 as fn depth
  bN <<- 10.85                 # coeff. for N0 as fn depth
  
} else if (flag_stn==2) {          # Biotrans (47 N 20 W)
  for (i in seq(1,13)) {
  MLDi[i] <- forcing[i,3]
  SSTi[i] <- forcing[i,4]
  }
  latitude <<- 47.0
  clouds <<- 6.0
  e0 <<- 12.0
  aN <<- 0.0174
  bN <<- 4.0
  
} else if (flag_stn==3) {          # Kerfix (50 40’S 68 25’E)
  for (i in seq(1,13)) {
  MLDi[i] <- forcing[i,7]
  SSTi[i] <- forcing[i,8]
  }
  latitude <<- -50.67
  clouds <<- 6.0
  e0 <<- 12.0
  aN <<- 0
  bN <<- 26.1
  
} else if (flag_stn==4) {          # Papa (50N 145W)
  for (i in seq(1,13)) {
  MLDi[i] <- forcing[i,5]
  SSTi[i] <- forcing[i,6]
  }
  latitude <<- 50.0
  clouds <<- 6.0
  e0 <<- 12.0
  aN <<- 0.0
  bN <<- 14.6
}

# Convert latitude to radians
latradians <<- latitude*pi/180.0

# ------------------------------------------------------------------------------- #
# Mixed layer depth: set up 366 unit array for one year                           #
# Interpolation used here from monthly data (assumed to be middle of month)       #
# ------------------------------------------------------------------------------- #

MLD <- array(dim=c(366))     # mixed layer depth (m)
SST <- array(dim=c(366))

for (istep in seq(1,365)) {   # interpolate data from monthly to daily values
rdiv <- 365.0/12.0            # monthly data
quot <- istep%/%rdiv          # quotient (integer)
remdr <- (istep%%rdiv)/rdiv   # remainder
n1 <- quot+1
if (n1==13) {
  n1 <- 1
}
n2 <- n1+1
iday <- istep+15          # set MLD data for middle (day 15) of month
if (iday>365) {
  iday <- iday-365
}
MLD[iday] <- MLDi[n1] + (MLDi[n2]-MLDi[n1])*remdr    # linear interpolation of data
SST[iday] <- SSTi[n1] + (SSTi[n2]-SSTi[n1])*remdr
}       # end loop: istep

for (istep in seq(366,2,by=-1)) {      # shift array because first position corresponds to t=0
MLD[istep] <- MLD[istep-1]
SST[istep] <- SST[istep-1]
}
MLD[1] <- MLD[366]           # note: position 1 in the array corresponds to t=0
SST[1] <- SST[366]
                            
######################################################
# -------------------------------------------------- #
# Variables specific to model: adjust accordingly    #
# -------------------------------------------------- #
######################################################

nSvar <- 4      # no. state variables
nfluxmax <- 10  # max no. of flux terms associated with any one state variable
nDvar <<- 20    # max no. auxilliary variables or output/plotting (these are defined in the Y array: FNget_flux function)

X <- c(Pinit,Ninit,Zinit,Dsinit)    # array X holds values of state variables: assign initial conditions

Svarname <- c("P","N","Z","D")      # array of text strings (keep them short) for state variables
Svarnames <- paste(Svarname[1],Svarname[2],Svarname[3],Svarname[4],sep=" ")   # concatenated string of state variable names


### End of variables specific to model  
 
#########################################################################
# --------------------------------------------------------------------- #
# Code below not to be changed for new models, except graph plotting    #
# --------------------------------------------------------------------- #
#########################################################################                       

# -------------------------- #
# Basic setup                #
# -------------------------- #

Y <<- rep(0.0, times=nDvar)                  # array Y holds values of auxiliary variables

ndays <- nyears*365                          # no. of days in simulation

nsteps <- floor(ndays/delta_t+delta_t/1000)                   # no. of time steps in simulation
nstepsday <- floor(1/delta_t+delta_t/1000)                    # no. of time steps per day
nwrite <- flag_outfreq*(ndays+1)+(1-flag_outfreq)*(nsteps+1)  # no. of times required for output to files and graphs
    
tm <- array(dim=c(nwrite), rep(0.0, times=nwrite))    # array to store time for output
flux <- matrix(0,nrow=nfluxmax,ncol=nSvar)            # matrix to hold update terms for state variables
fluxyear <- matrix(0,nrow=nfluxmax,ncol=nSvar)        # matrix to sum fluxes over year
Svar <- matrix(nrow=nwrite, ncol=nSvar)               # matrix to store values of state variables
Dvar <- matrix(nrow=nwrite, ncol=nDvar)               # matrix to store values of auxiliary variables

Svar[1,] <- X                                    # store initial values of state variables (position 1 in array corresponds to t=0)

tme <- 0       # time=0: start of simultion
tm[1] <- tme
 
Xout <- signif(X,digits=4)     # truncate X to 4 significant figures, ready for writing to output file
if (flag_outtype != 0) {
write(c(tme,Xout),file="out_statevars.txt",ncolumns=nSvar+1,sep=",",append=FALSE)  # write inital values of state vars to file
}
        
icount <- 1      # records position in array for writing output
istep <- 0

##############################    
# -------------------------- #
# Time loop starts here      #
# -------------------------- #
##############################

for (iyear in seq(1,nyears)) {             # loop over years

fluxyear <- rep(0.0,times=nfluxmax*nSvar)  # matrix to sum up fluxes over year

for (iday in seq(1,365)) {                 # loop over days of year

daynow <<- iday                            # global variable to pass to functions

for (istepnow in seq(1,nstepsday)) {       # loop over time steps within one day

istep <- istep + 1
istepplus <- istep+1               # (Dvar (auxiliary variables) stored at beginning of each time step)

fluxstep <- rep(0.0,times=nfluxmax*nSvar)  # matrix to sum up fluxes over time step
Ystep <- rep(0.0,times=nDvar)              # matrix to store average value of Y for each time step

# -------------------------- #
# Numerical integration      #
# -------------------------- # 

flux1 <- FNget_flux(X)*delta_t        # 1st iteration
deltaX1 <- colSums(flux1)             # sum of flux terms, for each state variable

if (flag_integ==0) {                  # Euler
  fluxstep <- fluxstep + flux1
  Ystep <- Ystep + Y
} else {                              # Runge Kutta 4
  fluxstep <- fluxstep + flux1/6
  Ystep <- Ystep + Y/6
}

if (flag_integ==1) {                  # Runge Kutta only

  flux2 <- FNget_flux(X + deltaX1/2)*delta_t       # 2nd iteration
  deltaX2 <- colSums(flux2)
  fluxstep <- fluxstep + flux2/3
  Ystep <- Ystep + Y/3
  
  flux3 <- FNget_flux(X + deltaX2/2)*delta_t       # 3rd iteration
  deltaX3 <- colSums(flux3)
  fluxstep <- fluxstep + flux3/3
  Ystep <- Ystep + Y/3
  
  flux4 <- FNget_flux(X + deltaX3)*delta_t         # 4th iteration
  deltaX4 <- colSums(flux4)
  fluxstep <- fluxstep + flux4/6
  Ystep <- Ystep + Y/6
}

# sum iterations
if (flag_integ==0) {              # Euler
  deltaX <- deltaX1
} else {    #Euler method  
  deltaX = (deltaX1 + 2*deltaX2 + 2*deltaX3 + deltaX4)/6    # Runge Kutta
}
    
X <- X + deltaX    # update state varialbes

fluxyear <- fluxyear + fluxstep

tme <- tme + delta_t    # update time

# -------------------------- #
# Write to output files      #
# -------------------------- # 

# 1st time step awkward: state vars already output for t=0; output fluxes and auxiliary variables for t=0 as calculated for first time step
if (istep==1) {
  Dvar[icount,] <- Ystep
 if (flag_outtype==2 | (flag_outtype==1 & nyears==1)) {                  # only needed if whole simulation from t=0 is to be written to files
    Dvarwrite <- signif(Dvar[icount,],digits=4)        # truncated to 4 significant figures: auxiliary variables
    fluxwrite <- signif(fluxstep,digits=4)             # truncated to 4 significant figures: fluxes
    tmenow <- 0.0     # actually, it is time 1 time step, but for the purposes of continuity in output files, set to zero
    write(c(tmenow,Dvarwrite),file="out_aux.txt",ncolumns=nDvar+1,sep=",",append=FALSE)
        for (iw in seq(1,nSvar)) {
    if (iw==1) {
      fluxwritecat <- fluxwrite[,1]
    } else {
      fluxwritecat <- c(fluxwritecat,fluxwrite[,iw])
    }
  }
    write(c(tmenow,fluxwritecat),file="out_fluxes.txt",ncolumns=nSvar*nfluxmax+1,sep=",",append=FALSE)
  }
}

# all other time steps
if (istepnow==nstepsday | flag_outfreq==0) {    # output on last time step of day (daily), or if every time step selected
  icount <- icount+1                            # records position in output arrays
  tm[icount] <- tme                             # time
  Svar[icount,] <- X                            # state variables
  Dvar[icount,] <- Ystep                        # auxiliary variables
  if ((flag_outtype==1 & (iyear==nyears | (iyear==nyears-1 & iday==365))) | flag_outtype==2) {
  appendtrue <- !(flag_outtype==1 & (iyear==nyears-1 & iday==365))   # wipe output file clear if output is for last year only, on last day of previous year
  Svarwrite <- signif(X,digits=4)               # truncated to  4 significant figures: state variables
  Dvarwrite <- signif(Dvar[icount,],digits=4)   # truncated to  4 significant figures: auxiliary variables     
  fluxwrite <- signif(fluxstep,digits=4)        # truncated to  4 significant figures: fluxes
  write(c(tme,Svarwrite),file="out_statevars.txt",ncolumns=nSvar+1,sep=",",append=appendtrue)
  write(c(tme,Dvarwrite),file="out_aux.txt",ncolumns=nDvar+1,sep=",",append=appendtrue)
  for (iw in seq(1,nSvar)) {
    if (iw==1) {
      fluxwritecat <- fluxwrite[,1]
    } else {
      fluxwritecat <- c(fluxwritecat,fluxwrite[,iw])
    }
  }
  write(c(tme,fluxwritecat),file="out_fluxes.txt",ncolumns=nSvar*nfluxmax+1,sep=",",append=appendtrue)
  }
}

# -------------------------- #
# Close loops: time          #
# -------------------------- # 

} # End loop: time steps over day
} # End loop: days of year    
} # End loop: years

# -------------------------------------------------------- #
# Print summed annual fluxes over last year to screen      #
# -------------------------------------------------------- # 

print("**************************************************")
hdr <- paste("annual fluxes (last year):",Svarnames,sep=" ")
print(hdr)
print("                                                  ")
print(fluxyear)

fluxvar <- colSums(fluxyear)    # fluxes for each state var summed over year
print("Total, each state var:                            ")
print(fluxvar)

print("                                                  ")

if (flag_outtype>0) {
  print("output: state varibles written to file out_statevars.txt")
  print("output: terms in differential equations written to file out_fluxes.txt")
  print("output: auxiliary variables written to file out_aux.txt")
} else {
  print("no output to text files (flag switched off)")  
}
print("                                                  ")
print("**************************************************")

# ----------------------------- #
# Graphs: Examples provided     #
# ----------------------------- # 

par(mfrow = c(3,3))    # layout for plots: 3 rows and 3 columns

# P, Z, D
maxy <- max(max(Svar[,1]),max(Svar[,3]),max(Svar[,4]))*1.2   # derive scaling for y on graph   
plot(c(0,ndays), c(0,maxy), xlab="days", ylab="mmol N m-3", type="n")
legend("topleft",legend=c("P","Z","D"),bty="n",lty=1,col=c("green","red","brown"),cex=1.2)
title(main="P, Z, D", col.main="red", font.main=4)
lines(tm,Svar[,1],col="green")
lines(tm,Svar[,3],col="red")
lines(tm,Svar[,4],col="brown")

# N
maxy <- max(Svar[,2])*1.2     
plot(c(0,ndays), c(0,maxy), xlab="days", ylab="mmol N m-3", type="n")
legend("topleft",legend=c("N"),bty="n",lty=1,col=c("blue"),cex=1.2)
title(main="N", col.main="red", font.main=4)
lines(tm,Svar[,2],col="blue")

# u (phytoplankton specific growth rate)
maxy <- max(Dvar[,3])*1.2     
plot(c(0,ndays), c(0,maxy), xlab="days", ylab="d-1", type="n")
legend("topleft",legend=c("u"),bty="n",lty=1,col=c("green"),cex=1.2)
title(main="u_P", col.main="red", font.main=4)
lines(tm,Dvar[,3],col="green")

# L_I, L_N
maxy <- max(max(Dvar[,1]),max(Dvar[,2]))   # derive scaling for y on graph   
plot(c(0,ndays), c(0,maxy), xlab="days", ylab="dimensionless", type="n")
legend("topleft",legend=c("L_I","L_N"),bty="n",lty=1,col=c("green","red"),cex=1.2)
title(main="L_I, L_N", col.main="red", font.main=4)
lines(tm,Dvar[,1],col="green")
lines(tm,Dvar[,2],col="red")

# MLD
maxy <- max(Dvar[,4])*1.2     
plot(c(0,ndays), c(0,maxy), xlab="days", ylab="m", type="n")
legend("topleft",legend=c("MLD"),bty="n",lty=1,col=c("blue"),cex=1.2)
title(main="MLD", col.main="red", font.main=4)
lines(tm,Dvar[,4],col="blue")

# noonPAR
maxy <- max(Dvar[,5])*1.2     
plot(c(0,ndays), c(0,maxy), xlab="days", ylab="W m-2", type="n")
legend("topleft",legend=c("noonPAR"),bty="n",lty=1,col=c("red"),cex=1.2)
title(main="noonpar", col.main="red", font.main=4)
lines(tm,Dvar[,5],col="red")

# chl
maxy <- max(Dvar[,6])*1.2     
plot(c(0,ndays), c(0,maxy), xlab="days", ylab="mg  m-3", type="n")
legend("topleft",legend=c("chl"),bty="n",lty=1,col=c("green"),cex=1.2)
title(main="chl", col.main="red", font.main=4)
lines(tm,Dvar[,6],col="green")





