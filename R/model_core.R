
#' Run the EMPOWER model
#'
#' @return Nothing (yet!)
#' @export
#'
EMPOWER <- function() {

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
  flag_stn    <<- init[9,1]     # choice of station; 1=India (60N 20W), 2=Biotrans (47N 20W), 3=Kerfix (50 40?S 68 25?E), 4=Papa (50N, 145W)
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

  } else if (flag_stn==3) {          # Kerfix (50 40?S 68 25?E)
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
}
