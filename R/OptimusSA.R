################################################################################
# Copyright (C) 2015+ Sahakyan Laboratory (aleksahak[at]cantab.net).           #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the following link:            #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################
#                                                                              #
#        Monte Carlo optimiser with Acceptance Ratio Annealing (MCARA)         #
#        Featuring:                                                            #
#        * an automatic adaptive pseudo-temperature bath control,              #
#        * memory friendly processing capability,                              #
#        * progress visualisation,                                             #
#        * multi-core = multi-replica optimisation.                            #
#                                                                              #
#*******************************************************************************
#-- NUMITER=6000     # Number of model optimisation steps.
#-- CYCLES = 2       # Number of annealing cycles (NUMITER should be divisable 
#                                                               to that number).
#*******************************************************************************
#-- ACCRATIO.IN=85   # Initial acceptance ratio (%) at the beginning of each
#                                                  constituent annealing cycles.
#-- ACCRATIO.FIN=5   # Final acceptance ratio (%) at the end of each constituent
#                                                              annealing cycles.
#-- STATWINDOW=70    # Number of last ongoing iterations to calculate acceptance
#                                         ratio for temperature auto-adjustment.
#*******************************************************************************
#-- T.INI=7.5        # Initial temperature (K) for Metropolis criterion.
#-- T.ADJSTEP=0.3    # Temperature change step-size for temperature auto-
#                               adjustment based on the actual acceptance ratio.
#-- TSCLnum=4        # Cutoff for one of the NumofAccRatSMIdeal and
#                          NumofAccRatGRIdeal numbers after which the adjustment
#                                               step is multiplied by T.SCALING.
#-- T.SCALING        # See above.
#-- T.MIN=0.001      # The minimum allowed temperature, set when the adjustments
#                                                               make T negative.
#*******************************************************************************
#-- DUMP.FREQ=10     # The frequency (in steps) of writing the found ongoing
#                                                                    best model.
#*******************************************************************************
#-- LIVEPLOT=TRUE    # Plotting the optimisation process in a pdf file.
#-- LIVEPLOT.FREQ=10 # Frequency (in steps) of plotting the results.
#-- PDFheight= 29    # Plot height in inches.
#-- PDFwidth = 20    # Plot width in inches.
#*******************************************************************************
#-- NCPU            # Number of CPU cores to use, by running more replicas of 
#                     the optimisation. The usage of more than 1 cores will
#                     attempt to load the foreach and doMC libraries.
#-- LONG = FALSE    # If TRUE, it means that a long simulation is expected to
#                     be done, hence the memory-friendly mode will be activated.
#-- SEED            # Setting the seed for the random number generator.
#-- OPTNAME         # The name of the optimisation process.
#*******************************************************************************
#--  Optimus() embeds a model to its core via 2 interfacing lines: a) by sourcing
#    MODEL.PATH (local=FALSE) R file at the initialisation part of mcara(), and
#    b) by sourcing (local=TRUE) the MDLCALL.PATH R file from inside the optimi-
#    sation loop, where the MDLCALL.PATH file is supposed to use the definitions
#    specified in MODEL.PATH to generate an instance of the mdl object. The
#    latter consists of the components $sol (a vector of prediction instance),
#    $coef.new (a vector of coefficient instance) and $sol.perf (an object of
#    the solution performance with $Q /quality score, such as RMSE etc,/ and $E
#    /pseudo-energy/ components).
#-- MODEL.PATH      # Path to the required R file that holds the model details.
#*******************************************************************************
################################################################################
Optimus <- function(NUMITER       = 1000000,
                    CYCLES        = 10,
                    ACCRATIO.IN   = 90,
                    ACCRATIO.FIN  = 0.5,
                    STATWINDOW    = 70,
                    T.INI         = 0.00001,
                    T.ADJSTEP     = 0.000000005,
                    TSCLnum       = 2,
                    T.SCALING     = 3,
                    T.MIN         = 0.000000005,
                    DUMP.FREQ     = 10000,
                    LIVEPLOT      = TRUE,
                    LIVEPLOT.FREQ = 100000,
                    PDFheight     = 29,
                    PDFwidth      = 20,
                    NCPU          = 4,  
                    LONG          = TRUE,
                    SEED          = 840,
                    OPTNAME       = "poly",
                    MODEL.PATH    = "OptMDL_1.R"
                   ){
################################################################################

  set.seed(SEED)
  suppressWarnings(suppressPackageStartupMessages(library(foreach)))
  if(NCPU > 1){
    suppressWarnings(suppressPackageStartupMessages(library(doMC)))
    registerDoMC(cores = NCPU)
    `%op%` <- `%dopar%`
    print(paste("Running ", NCPU, " replicas of optimisation.", sep=""), quote=F)
  } else {
    `%op%` <- `%do%`
  }
  
  seeds <- floor(runif(n=NCPU, min=1, max=1000))

  #-- INITIALISATION # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                              #-- The target energies of all the NUMITER trials.
  ENERGY.TRIAL.VEC      = NULL 
             #-- The energy of the system, updated only if the move is accepted.
  ENERGY.DE.FACTO       = NULL 
                                    #-- Probability of acceptance for each step.
  PROB.VEC              = NULL 
                 #-- Acceptance vector for each step: 1/0 = accept/don't accept.
  ACCEPTANCE            = NULL 
                                         #-- The real temperature for each step.
  T.DE.FACTO            = NULL 
                                           #-- The aimed ideal acceptance ratio.
  IDEAL.ACC.VEC         = rep(seq(from = ACCRATIO.IN, to = ACCRATIO.FIN, by = 
                                 (ACCRATIO.FIN-ACCRATIO.IN)/((NUMITER/CYCLES)-1)
                                 ), CYCLES) 
           #-- The calculated de facto acceptance rate over the last STATWINDOW.
  ACC.VEC.DE.FACTO      = NULL 
                                         #-- The corresponding checkpoint steps.
  STEP4ACC.VEC.DE.FACTO = NULL 
                                              #-- The quality of the prediction.
  Q.STRG                = NULL
                                                            #-- The stored step.
  STEP.STORED           = NULL

          #-- initial big value, so that the first step will always be accepted!
  E.old         <- 100000  
                                  #-- initial T for starting further adjustments 
  T             <- T.INI    
                                                  #-- initial meaningless value!
  Q.old         <- 0
                        #-- initial big value for the lowest energy model stored        
  E.stored      <- 100000
  Step.stored   <- 1
  STEP.add      <- 1

         #-- Number of times acceptance ratio is sequentially smaller than ideal
  NumofAccRatSMIdeal <- 0
         #-- Number of times acceptance ratio is sequentially greater than ideal
  NumofAccRatGRIdeal <- 0
  t.adjstep        <- T.ADJSTEP
  AccR.category    <- "INITIAL"
  new.T.INI        <- T.INI
  instanceOFswitch <- 0
  #^^ INITIALISATION # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



  print("Monte Carlo optimisation with acceptance annealing...", quote=FALSE)
  
 #-- PARALLEL PROCESSING WRAP # # # # # # # # #
 foreach(repl=1:NCPU, .inorder=FALSE) %op% {
 
  set.seed(seeds[repl])
  if(NCPU==1){ repl <- NULL }
  
  ###***************
  ### MODEL ########
  source(MODEL.PATH)
  ##^ MODEL ########
  ###***************
  
  #-- STARTING MC ITERATIONS # # # # # # # # # # # # # # # # # # # # # # # # # #
  for(STEP in 1:NUMITER){  #-- step of the trial, from 1 to NUMITER

    if(STEP%%1000==0){
      print(paste("Step ",STEP," out of ",NUMITER,sep=""), quote=FALSE)
    }
    
    ###***************
    ### MODEL ########
    K.new <- r(K=K)
    O     <- m(K=K.new)
    snp   <- u(O=O)
    ##^ MODEL ########
    ###***************    
  
    Q       <- snp$Q
    E       <- snp$E                      ; ENERGY.TRIAL.VEC[STEP.add]<- E
    DE      <- E - E.old                  ; T.DE.FACTO[STEP.add]      <- T
    #-- Metropolis criterion
    P.accept<- min(1, exp(-DE/T))         ; PROB.VEC[STEP.add]        <- P.accept
    accept  <- rbinom(1, 1, prob=P.accept); ACCEPTANCE[STEP.add]      <- accept

    ################
    if( accept==1 ){

      E.old  <- E
      K      <- K.new
      Q.old  <- Q

      #-- # # # # # # # # #
      if(E.old < E.stored){
        E.stored       <- E.old
        Step.stored    <- STEP
        K.stored       <- K
        O.stored       <- O
        DUMP.MODEL     <- NULL
        DUMP.MODEL[1]  <- "QUALITY:"
        DUMP.MODEL[2]  <- paste("E: ",round(E.old,3),sep="")
        DUMP.MODEL[3]  <- paste("Q: ",round(Q.old,3),sep="")
#        DUMP.MODEL[4]  <- "TERMS:"
#        DUMP.MODEL[5]  <- paste(as.character(names(K.stored)), collapse=" ")
#        DUMP.MODEL[6]  <- "COEFFICIENTS:"
#        DUMP.MODEL[7]  <- paste(format(as.vector(K.stored), scientific=F, trim=T), collapse=" ")
#        DUMP.MODEL[8]  <- "OBSERVABLES:"
#        DUMP.MODEL[9]  <- paste(as.character(names(O.stored)), collapse=" ")
#        DUMP.MODEL[10] <- "PREDICTIONS:"
#        DUMP.MODEL[11] <- paste(format(as.vector(O.stored), scientific=F, trim=T), collapse=" ")
#        DUMP.MODEL[12] <- "TARGET:"
#        DUMP.MODEL[13] <- paste(format(as.vector(y), scientific=F, trim=T), collapse=" ")
#        DUMP.MODEL[14] <- "################################"
      }
      #-- # # # # # # # # #

    }
    ################

    ENERGY.DE.FACTO[STEP.add] <- E.old
    Q.STRG[STEP.add]          <- Q.old
    STEP.STORED[STEP.add]     <- STEP

    # NumofAccRatSMIdeal -- Number of times acceptance ratio is sequentially smaller than ideal
    # NumofAccRatGRIdeal -- Number of times acceptance ratio is sequentially greater than ideal
    # TSCLnum=4          -- Cutoff for one of the above numbers after which the adjustment step
    #                       for the temperature will be twice increased or decreased.
    # t.adjstep <- T.ADJSTEP
    # AccR.category

    #-- Adjusting the temperature every STATWINDOW STEP
    if((length(ACCEPTANCE)%%STATWINDOW) == 0){ #-- thus the number of entries is n*STATWINDOW
      accept.ratio.inlastWIN <- 100*sum(ACCEPTANCE[(length(ACCEPTANCE)-STATWINDOW+1):length(ACCEPTANCE)])/
                                                          STATWINDOW
      ACC.VEC.DE.FACTO <- c(ACC.VEC.DE.FACTO, accept.ratio.inlastWIN)
      STEP4ACC.VEC.DE.FACTO <- c(STEP4ACC.VEC.DE.FACTO, STEP)

      if(accept.ratio.inlastWIN < IDEAL.ACC.VEC[STEP]){ #-- CATEGORY - SMALLER
        NumofAccRatSMIdeal <- NumofAccRatSMIdeal + 1
        NumofAccRatGRIdeal <- 0
        if(NumofAccRatSMIdeal==TSCLnum){ 
          t.adjstep <- T.SCALING*t.adjstep; NumofAccRatSMIdeal <- 0
        }
        if(AccR.category=="GREATER"){
          t.adjstep <- T.ADJSTEP; instanceOFswitch <- instanceOFswitch + 1 
          if(instanceOFswitch==1){ new.T.INI <- T } 
        } # resets t.adjstep, if it has just crossed the ideality line.
        T <- T + t.adjstep
        AccR.category <- "SMALLER"
      }
      if(accept.ratio.inlastWIN >= IDEAL.ACC.VEC[STEP]){ #-- CATEGORY - GREATER
        NumofAccRatSMIdeal <- 0
        NumofAccRatGRIdeal <- NumofAccRatGRIdeal + 1
        if(NumofAccRatGRIdeal==TSCLnum){
          t.adjstep <- T.SCALING*t.adjstep; NumofAccRatGRIdeal <- 0
        }
        if(AccR.category=="SMALLER"){
          t.adjstep <- T.ADJSTEP; instanceOFswitch <- instanceOFswitch + 1
          if(instanceOFswitch==1){ new.T.INI <- T }
        } # resets t.adjstep, if it has just crossed the ideality line.
        T <- T - t.adjstep
        #-- not allowing minus temperatures (absolute scale) and 0 (DE/T=NaN)
        if(T <= 0){ T <- T.MIN }
        AccR.category <- "GREATER"
      }
    }

    #-- reseting temperature after a cycle is completed
    if(STEP%%(NUMITER/CYCLES)==0){
        T <- new.T.INI; instanceOFswitch <- 0
    }
    #-- # # # # # # # # # # # # # # # # # # # # # # # #

    #-- PLOTTING  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    if(LIVEPLOT==TRUE){
      if(STEP%%LIVEPLOT.FREQ == 0){
      
        pdf(width=PDFwidth, height=PDFheight, file=paste(OPTNAME,repl,".pdf",sep=""))
        par(mfrow=c(5,1))

        plot(y=PROB.VEC, x=STEP.STORED, ylab="Acceptance P", xlab="Step", 
             ylim=c(0,1), xlim=range(STEP.STORED), type="l", col="blue", lwd=1)
        lines(y=c(0,1), x=rep(Step.stored,2), lwd="2.5", col="green")

        plot(y=T.DE.FACTO, x=STEP.STORED, ylab="Temperature", xlab="Step",
             xlim=range(STEP.STORED), type="l", col="orange", lwd=2)
        lines(y=range(T.DE.FACTO), x=rep(Step.stored,2), lwd="2.5", col="green")

        plot(y=IDEAL.ACC.VEC[STEP.STORED], x=STEP.STORED,
             ylab="Acceptance ratios (%)", xlab="Step", ylim=c(0,100),
             xlim=range(STEP.STORED), type="l",lty="dashed", lwd=2)
        lines(y=c(0,100), x=rep(Step.stored,2), lwd="2.5", col="green")
        if(length(ACC.VEC.DE.FACTO)!=0){
          lines(y=ACC.VEC.DE.FACTO, x=STEP4ACC.VEC.DE.FACTO, ylim=c(0,100),
                xlim=range(STEP.STORED), col="red", lwd=2)
        }

        plot(y=ENERGY.DE.FACTO, x=STEP.STORED, ylab="Energy", xlab="Step",
             xlim=range(STEP.STORED), type="l", col="navy", lwd=2)
        lines(y=range(ENERGY.DE.FACTO), x=rep(Step.stored,2), lwd="2.5",
              col="green")

        plot(y=Q.STRG, x=STEP.STORED, ylab="Q", xlab="Step",
             xlim=range(STEP.STORED), type="l", col="blue", lwd=2)
        lines(y=range(Q.STRG), x=rep(Step.stored,2), lwd="2.5", col="green")

        dev.off()
        
      }
    }
    #^^ PLOTTING  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


    if( STEP%%DUMP.FREQ == 0 & exists("DUMP.MODEL") ){
      write(DUMP.MODEL, file=paste(OPTNAME,repl,"_model_QE.log",sep=""))
      save(K.stored,    file=paste(OPTNAME,repl,"_model_K.Rdata",sep=""))
      save(O.stored,    file=paste(OPTNAME,repl,"_model_O.Rdata",sep=""))
    }

    STEP.add <- STEP.add + 1    ########

    #-- Trimming the data holding vectors to be memory-friendly if LONG==TRUE
    if(LONG==TRUE){
      if(STEP%%10000 == 0 & STEP!=10000){
        ENERGY.DE.FACTO       <- ENERGY.DE.FACTO[(length(ENERGY.DE.FACTO)-10000+1):length(ENERGY.DE.FACTO)]
        Q.STRG                <- Q.STRG[(length(Q.STRG)-10000+1):length(Q.STRG)]
        STEP.STORED           <- STEP.STORED[(length(STEP.STORED)-10000+1):length(STEP.STORED)]
        ENERGY.TRIAL.VEC      <- ENERGY.TRIAL.VEC[(length(ENERGY.TRIAL.VEC)-10000+1):length(ENERGY.TRIAL.VEC)]
        T.DE.FACTO            <- T.DE.FACTO[(length(T.DE.FACTO)-10000+1):length(T.DE.FACTO)]
        PROB.VEC              <- PROB.VEC[(length(PROB.VEC)-10000+1):length(PROB.VEC)]
        ACCEPTANCE            <- ACCEPTANCE[(length(ACCEPTANCE)-10000+1):length(ACCEPTANCE)]
        ACC.VEC.DE.FACTO      <- ACC.VEC.DE.FACTO[(length(ACC.VEC.DE.FACTO)-(10000/STATWINDOW)+1):length(ACC.VEC.DE.FACTO)]
        STEP4ACC.VEC.DE.FACTO <- STEP4ACC.VEC.DE.FACTO[(length(STEP4ACC.VEC.DE.FACTO)-(10000/STATWINDOW)+1):length(STEP4ACC.VEC.DE.FACTO)]
        STEP.add              <- 10001
      }
    }
    #^^ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  } ##-- for(STEP in 1:NUMITER)
  #^^ STARTING MC ITERATIONS # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  #-- if LONG==TRUE, OUTPUT will only hold the trimmed data!
  OUTPUT                       <- NULL
  OUTPUT$K.stored              <- K.stored
  OUTPUT$O.stored              <- O.stored
  OUTPUT$STEP                  <- STEP
  OUTPUT$PROB.VEC              <- PROB.VEC
  OUTPUT$T.DE.FACTO            <- T.DE.FACTO
  OUTPUT$IDEAL.ACC.VEC         <- IDEAL.ACC.VEC
  OUTPUT$ACC.VEC.DE.FACTO      <- ACC.VEC.DE.FACTO
  OUTPUT$STEP4ACC.VEC.DE.FACTO <- STEP4ACC.VEC.DE.FACTO
  OUTPUT$ENERGY.DE.FACTO       <- ENERGY.DE.FACTO
  OUTPUT$Q.STRG                <- Q.STRG
  OUTPUT$ACCEPTANCE            <- ACCEPTANCE
  OUTPUT$STEP.STORED           <- STEP.STORED

  save(OUTPUT, file=paste(OPTNAME,repl,"_model_ALL.Rdata", sep=""))
 
 }
 #^^ PARALLEL PROCESSING WRAP # # # # # # # # #
  
  print("An optimal model is obtained and the data are saved !!!", quote=FALSE)
  if(NCPU==1){ return(OUTPUT) } else { return(0) }
  
}
################################################################################
