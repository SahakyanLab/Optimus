################################################################################
# Shuffling a set of genomic contacts to form a new set using Optimus 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
wk.dir = "/Optimus/"
data.dir = paste0(wk.dir, "/documentation/example_5/TutorialData")
out.dir = paste0(wk.dir, "/documentation/example_5")
setwd(out.dir)
### OTHER SETTINGS #############################################################
region.len = 4e+04
mingap.Mb = 2e+06

# Optimus SA | RE
opt.type = "SA" # "SA" | # RE
# If opt.type = "RE"
ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2)

out.name = "IJ.NEW.OPTI.SA" # "IJ.NEW.OPTI.SA" | "IJ.NEW.OPTI.RE"
replica = 4 # 4 | # 12

numiter = 2e+05 
cycles = 2
dump.freq = 1e+05
calcErrorOnly = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
source(paste0(wk.dir, "/R/Optimus.R"))
source(paste0(wk.dir, "/R/OptimusSA.R"))
source(paste0(wk.dir, "/R/TemperatureControlUnit.R"))

# Function to calculate error rate of shuffled set of contacts in comparison to
# original set and based on the following criteria for a valid contact:
# 1. i should always be greater than j.  
# 2. i and j should be separated by a distance greater than 2 Mb. 
# 3. The set should only contain unique pairs/contacts.
# 4. Contact in new set should not be in the original set. 
calcErrorRate <- function(
  test = IJ.NEW,
  ref = IJ.ORIG,
  # Should be > this value; in terms of bins
  gap = 50
){
  
  # Not satisfying gap | present in orig set | duplicated in new set
  x <- sum( 
    ( (test[,"j"]-test[,"i"]) <= gap )  | (test[,"j"]==ref[,"j"])  | duplicated(test) 
  )
  
  return(x/nrow(IJ.ORIG)*100)
  
}

# Define interfacing functions for Optimus
u <- function(O, DATA = NULL){
  result <- NULL
  result$Q <- -O
  result$E <- O
  return(result)
}

r <- function(K){
  K.new <- K
  new.ind <- sample(x=1:length(K.new[,"j"]), size=2, replace=FALSE)
  K.new[new.ind,"j"] <- K.new[rev(new.ind),"j"]
  return(K.new)
}

m <- function(K, DATA = NULL){
  errorCont <- sum( ( (K[,"j"]-K[,"i"]) <= DATA$gaplimit )  | ( K[,"j"]==DATA$IJ.ORIG[,"j"] )  | duplicated(K) )
  O <- (errorCont/DATA$numContacts)*100
  return(O)
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
start.time <- Sys.time() 

# Required gap between contacting bins; in terms of bins
mingap <- mingap.Mb/region.len

# Load IJ.ORIG
load(file=paste0(data.dir, "/ij_orig.RData"))

if(calcErrorOnly==FALSE){
  
  DATA <- list(IJ.ORIG=IJ.ORIG, numContacts=nrow(IJ.ORIG), gaplimit=mingap)
  
  K <- IJ.ORIG
  # First random shuffle of js; is are kept in the same order
  K[,"j"] <- sample(x=K[,"j"], size=nrow(K), replace=FALSE)
  
  # Uses default SEED = 840
  
  if(opt.type=="SA"){
    Optimus(NCPU = replica, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = opt.type,
            OPTNAME = out.name, DATA = DATA, NUMITER = numiter, CYCLES = cycles, DUMP.FREQ = dump.freq,
            LONG = TRUE)
  } else if(opt.type=="SA"){
    # Uses default SEED = 840
    Optimus(NCPU = replica, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = ACCRATIO,
            OPT.TYPE = opt.type, DATA = DATA, OPTNAME = out.name, NUMITER = numiter, STATWINDOW = 50,
            DUMP.FREQ = dump.freq, LONG = TRUE)
  } else {
    stop("opt.type can only be either SA or RE")
  }
  
  rm(DATA, K);gc()
    
} 

#--------------------------------------
# Calculate error rate for all replicas
#---------------------------------------
ePERC.v <- rep(NA, times=replica)

for(x in 1:replica){
  
  if(opt.type=="SA"){
    load(file=paste0(out.dir, "/", out.name, x, "_model_K.Rdata"))
    test.df <- K.stored
  } else if(opt.type=="RE"){
    load(file=paste0(out.dir, "/", out.name, x, "_model_ALL.Rdata"))
    test.df <- OUTPUT$K.stored
  }
    
  ePERC <- calcErrorRate(
    test=test.df,
    ref=IJ.ORIG,
    # Should be > this value; in terms of bins
    gap=mingap
  )
  ePERC.v[x] <- paste0("%ERROR", x, ": ", round(x=ePERC, digits=10))
}
write(ePERC.v, file=paste0(out.dir, "/stats_", out.name))

end.time <- Sys.time()
end.time-start.time

# rm(list=ls())
