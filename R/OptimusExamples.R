#' @export

OptimusExamples <- function(example=1,dir='.',file_name='example.R',method='SA',mopac="~/Downloads/MOPAC2016_for_Macintosh/MOPAC2016.exe",run=FALSE){
  setwd(dir)
  dir.create(paste0('example_',example))
  setwd(paste0(dir,'/example_',example))
  getwd()
  if (example==1) {
    text = 'library(Optimus)
set.seed(845)
x <- runif(1000, min=-15, max=10)
y <- -1.0*x - 0.3*x^2 + 0.2*x^3 + 0.01*x^4 + rnorm(length(x),mean=0,sd=30)

DATA   <- NULL
DATA$x <- x
DATA$y <- y

################################################################################
K <- c(k1=1.0, k2=1.0, k3=1.0, k4=1.0)
################################################################################

################################################################################
m <- function(K,DATA){
  x <- DATA$x
  O <- K["k1"]*x + K["k2"]*x^2 + K["k3"]*x^3 + K["k4"]*x^4
  return(O)
}
################################################################################

################################################################################
u <- function(O,DATA){
  y <- DATA$y
  Q <- sqrt(mean((O-y)^2))
  E <- Q # For RMSD, <-> negative sign or other mathematical operation
  # is not needed.
  
  RESULT <- NULL
  RESULT$Q <- Q
  RESULT$E <- E
  return(RESULT)
}
################################################################################

################################################################################
r <- function(K){
  K.new <- K
  # Randomly selecting a coefficient to alter:
  K.ind.toalter <- sample(size=1, x=1:length(K.new))
  # Creating a potentially new set of coefficients where one entry is altered
  # by either +move.step or -move.step, also randomly selected:
  move.step <- 0.0005
  K.new[K.ind.toalter] <- K.new[K.ind.toalter] + sample(size=1, x=c(-move.step, move.step))
  
  ## Setting the negative coefficients to 0 (not necessary in this example,
  ## useful for optimising rate constants):
  #neg.ind <- which(K.new < 0)
  #if(length(neg.ind)>0){ K.new[neg.ind] <- 0 }
  
  return(K.new)
}
################################################################################
    '
    if (method=='SA') {
      call = 
        'Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", DATA = DATA, OPTNAME = "poly_4_SA", LONG = FALSE)'
    } else if (method=='RE') {
      call = 
        'Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2), OPT.TYPE = "RE", DATA = DATA, OPTNAME = "poly_12_RE", LONG = FALSE)'
    }
    
    
  } else if (example==2) {
    text <- 'library(Optimus)
set.seed(845)
x <- runif(1000, min = -15, max = 10)
y <- -1 * x - 0.3 * x^2 + 0.2 * x^3 + 0.01 * x^4 + rnorm(length(x), mean = 0, sd = 30)
DATA <- NULL
DATA$x <- x
DATA$y <- y

################################################################################
K <- c(term1=rbinom(n=1, size=1, prob=0.5),
       term2=rbinom(n=1, size=1, prob=0.5),
       term3=rbinom(n=1, size=1, prob=0.5),
       term4=rbinom(n=1, size=1, prob=0.5),
       term5=rbinom(n=1, size=1, prob=0.5),
       term6=rbinom(n=1, size=1, prob=0.5),
       term7=rbinom(n=1, size=1, prob=0.5),
       term8=rbinom(n=1, size=1, prob=0.5),
       term9=rbinom(n=1, size=1, prob=0.5),
       term10=rbinom(n=1, size=1, prob=0.5),
       term11=rbinom(n=1, size=1, prob=0.5),
       term12=rbinom(n=1, size=1, prob=0.5),
       term13=rbinom(n=1, size=1, prob=0.5),
       term14=rbinom(n=1, size=1, prob=0.5),
       term15=rbinom(n=1, size=1, prob=0.5),
       term16=rbinom(n=1, size=1, prob=0.5),
       term17=rbinom(n=1, size=1, prob=0.5),
       term18=rbinom(n=1, size=1, prob=0.5),
       term19=rbinom(n=1, size=1, prob=0.5),
       term20=rbinom(n=1, size=1, prob=0.5),
       term21=rbinom(n=1, size=1, prob=0.5),
       term22=rbinom(n=1, size=1, prob=0.5),
       term23=rbinom(n=1, size=1, prob=0.5),
       term24=rbinom(n=1, size=1, prob=0.5),
       term25=rbinom(n=1, size=1, prob=0.5),
       term26=rbinom(n=1, size=1, prob=0.5),
       term27=rbinom(n=1, size=1, prob=0.5),
       term28=rbinom(n=1, size=1, prob=0.5),
       term29=rbinom(n=1, size=1, prob=0.5),
       term30=rbinom(n=1, size=1, prob=0.5))
################################################################################

################################################################################
m <- function(K, DATA){
  y <- DATA$y
  x <- DATA$x
  
  terms <- c("+x",
             "+I(x^2)",
             "+I(x^3)",
             "+I(x^4)",
             "+I(x^5)",
             "+I(x^6)",
             "+I(x^7)",
             "+I(x^8)",
             "+I(x^9)",
             "+I(x^10)",
             "+I(exp(x))",
             "+I(abs(x))",
             "+I(sin(x))",
             "+I(cos(x))",
             "+I(tan(x))",
             "+I(sin(x)*cos(x))",
             "+I((sin(x))^2)",
             "+I((cos(x))^2)",
             "+I(sin(x^2))",
             "+I(sin(x^3))",
             "+I(cos(x^2))",
             "+I(cos(x^3))",
             "+I(sin(x^3)*cos(-x))",
             "+I(cos(x^3)*sin(-x))",
             "+I(sin(x^5)*cos(-x))",
             "+I(cos(x^5)*sin(-x))",
             "+I(exp(x)*sin(x))",
             "+I(exp(x)*cos(x))",
             "+I(abs(x)*sin(x))",
             "+I(abs(x)*cos(x))")
  
  # Intercept term is allowed in all the cases of the equation.
  ind.terms <- which(K == 1)
  if(length(ind.terms)!=0){
    equation <- paste(c("y~",terms[ind.terms]), collapse="")
  } else {
    equation <-"y~x" # In case there are no active terms, use a simple linear model.
  }
  
  O <- glm(equation, data = environment())
  
  return(O)
}
################################################################################

################################################################################
u <- function(O, DATA){
  y <- DATA$y
  
  Q <- sqrt(mean((O$fitted.values-y)^2))
  E <- AIC(O)/1000 # Akaike\'s information criterion.
  
  result   <- NULL
  result$Q <- Q
  result$E <- E
  return(result)
}
################################################################################

################################################################################
r <- function(K){
  K.new <- K
  # Randomly selecting a term:
  K.ind.toalter <- sample(size=1, x=1:length(K.new))
  # If the term is on (1), switching it off (0) or vice versa:
  if(K.new[K.ind.toalter]==1){
    K.new[K.ind.toalter] <- 0
  } else {
    K.new[K.ind.toalter] <- 1
  }
  return(K.new)
}
################################################################################
'
    if (method=='SA') {
      call <- 
        'Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", DATA = DATA, OPTNAME = "term_4_SA", NUMITER = 200000, CYCLES = 2, DUMP.FREQ = 100000, LONG = FALSE)'
    } else if (method=='RE') {
      call <-
        'Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2), OPT.TYPE = "RE", DATA = DATA, OPTNAME = "term_12_RE", NUMITER = 200000, STATWINDOW = 50, DUMP.FREQ = 100000, LONG = FALSE)'
    }
    
    
  } else if (example==3) {
    text <- paste0('library(Optimus)
set.seed(845)
################################################################################
K <- c(PHI=90, PSI=90)
################################################################################

################################################################################
m <- function(K,
              notconvergedE = -100.00,
              mopac.cmd="',mopac,'",
              clean = TRUE){
  # mopac.cmd="~/Downloads/MOPAC2016_for_Macintosh/MOPAC2016.exe"
  # m(K=c(PHI=90.0, PSI=90.0), notconvergedE = -100.00,
  #   mopac.cmd="/home/alex/prog/mopac2016/MOPAC2016.exe", clean=TRUE)
  
  
  # MOPAC semiempirical QM input file preparation, with given PHI and PSI
  # dihedral angles.
  
  geo <- c(
    "RHF PM6 EF GEO-OK MMOK T=10 THREADS=1",
    "Vitamin C with two controllable dihedral angles psi(7,6,3,1) and phi(11,10,6,7)",
    "  ",
    "O     0.00000000  0    0.0000000  0    0.0000000  0     0     0     0",
    "H     0.98468620  1    0.0000000  0    0.0000000  0     1     0     0",
    "C     1.43651250  1  110.7230618  1    0.0000000  0     1     2     0",
    "H     1.10751723  1  103.6603154  1 -167.5282722  1     3     1     2",
    "H     1.10658657  1  110.2236860  1  -51.3620456  1     3     1     2",
    "C     1.53950336  1  112.8074046  1 -123.2791585  1     3     4     5",
    paste0("O     1.42824262  1  103.4315186  1 ",K["PSI"]," 0     6     3     1"),
    "H     0.99584949  1  109.9022382  1 -165.7055126  1     7     6     3",
    "H     1.11472171  1  108.4417082  1   75.1535637  1     6     7     8",
    "C     1.54244170  1  109.4042184  1 -120.8240216  1     6     7     9",
    paste0("O     1.46313669  1  105.7792445  1 ",K["PHI"]," 0    10     6     7"),
    "H     1.11252563  1  112.8336666  1 -114.5813834  1    10     6    11",
    "C     1.51686608  1  113.4849244  1 -112.8332453  1    10    12    11",
    "O     1.34410484  1  125.3617342  1  179.6090511  1    13    10    11",
    "H     1.03381724  1  110.9736522  1  -13.3419919  1    14    13    10",
    "C     1.36084908  1  124.8906459  1  167.6242325  1    13    14    15",
    "O     1.35614887  1  131.9374989  1   -0.0333000  1    16    13    14",
    "H     1.00338885  1  109.4220239  1    0.3798200  1    17    16    13",
    "C     1.49109250  1  118.0837177  1 -179.7749947  1    16    17    18",
    "O     1.18961787  1  136.9144035  1   -0.6060924  1    19    16    17",
    "  "
  )
  
  # Submitting the MOPAC optimisation job, where all the spatial parameters
  # are relaxed except the pre-set PHI and PSI angles. The job is run requesting
  # maximum 10 seconds of time limitation. Most (if not all) complete within
  # half a second.
  random.id <- as.character(sample(size=1, x=1:10000000))
  write(geo, file=paste0(random.id,".mop"))
  system(paste0(mopac.cmd," ",random.id,".mop"))
  
  if( file.exists(paste0(random.id,".arc")) ){
    e.line <- grep("HEAT OF FORMATION",
                   readLines(paste0(random.id,".arc")),
                   value=TRUE)
    e.line <-  strsplit(e.line," ")[[1]]
    O <- as.numeric(e.line[e.line!=""][5])
  } else {
    O <- notconvergedE
  }
  
  if(clean){
    file.remove(grep(random.id, dir(), value=TRUE))
  }
  
  return(O) # heat of form ation in kcal/mol
}
################################################################################

################################################################################
u <- function(O){
  result   <- NULL
  result$Q <- -O
  result$E <- O
  return(result)
}
################################################################################

################################################################################
r <- function(K){
  K.new <- K
  # Setting the alteration angle to 3 degrees:
  alter.by <- 3
  # Randomly selecting a term:
  K.ind.toalter <- sample(size=1, x=1:length(K.new))
  # Altering that term by either +alter.by or -alter.by
  K.new[K.ind.toalter] <-
    K.new[K.ind.toalter] + sample(size=1, x=c(alter.by, -alter.by))
  
  # Setting the dihedral angles to be always within the -180 to 180 range.
  if( K.new[K.ind.toalter] > 180.0 ){
    K.new[K.ind.toalter] <- K.new[K.ind.toalter] - 360
  }
  
  if( K.new[K.ind.toalter] < -180.0 ){
    K.new[K.ind.toalter] <- K.new[K.ind.toalter] + 360
  }
  
  return(K.new)
}
################################################################################
')
    if (method=='SA') {
      call <- 
      'Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", OPTNAME = "vitamin_4_SA", NUMITER = 1e+05, CYCLES = 2, DUMP.FREQ = 50000, LONG = FALSE)'
    } else if (method=='RE') {
      call <- 
      'Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2), OPT.TYPE = "RE", OPTNAME = "vitamin_12_RE", NUMITER = 1e+05, EXCHANGE.FREQ = 500, STATWINDOW = 50, DUMP.FREQ = 50000, LONG = FALSE)'
    }
    
  } else if (example==4) {
    text <- 'library(Optimus)
set.seed(845)
state  <- c(cA=100, cB=100, cC=100, cAP=0, cBP=0, cCP=0)
target <- c(cA=90, cB=20, cC=70, cAP=10, cBP=80, cCP=30)
model <- function(t, state, K){
  
  with( as.list(c(state, K)), {
    # rate of change
    dcA  <- -k1*cA+k2*cAP*cB
    dcB  <- -k2*cAP*cB+k3*cBP*cC
    dcC  <- -k3*cBP*cC+k4*cCP
    dcAP <- -dcA
    dcBP <- -dcB
    dcCP <- -dcC
    # return the rate of change
    list(c(dcA, dcB, dcC, dcAP, dcBP, dcCP))
  })
}

DATA <- NULL
DATA$state  <- state
DATA$target <- target
DATA$model  <- model

################################################################################
K <- c(k1=1.0, k2=1.0, k3=1.0, k4=1.0)
################################################################################

################################################################################
library(deSolve)
m <- function(K, DATA){
  state <- DATA$state
  model <- DATA$model
  
  span = 10.0
  
  times <- c(0, span)
  O   <- ode(y=state, times=times, func=model, parms=K)[2,2:(length(state)+1)]
  return(O)
}
################################################################################

################################################################################
u <- function(O, DATA){
  target <- DATA$target
  RESULT <- NULL
  RESULT$Q <- sqrt(mean((O-target)^2)) # measure of the fit quality
  RESULT$E <- RESULT$Q # the pseudoenergy derived from the above measure
  
  return(RESULT)
}
################################################################################

################################################################################
r <- function(K){                               
  K.new <- K
  # Randomly selecting a coefficient to alter:
  K.ind.toalter <- sample(size=1, x=1:length(K.new))
  # Creating a potentially new set of coefficients where one entry is altered
  # by either +move.step or -move.step, also randomly selected:
  move.step <- 0.0002
  K.new[K.ind.toalter] <- K.new[K.ind.toalter] + sample(size=1, x=c(-move.step, move.step))
  
  ## Setting the negative coefficients to 0
  neg.ind <- which(K.new < 0)
  if(length(neg.ind)>0){ K.new[neg.ind] <- 0 }
  
  return(K.new)
}
################################################################################
'
    if (method=='SA') {
      call <- 
      'Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", OPTNAME = "DE_4_SA", DATA = DATA, NUMITER = 200000, CYCLES = 2, DUMP.FREQ = 100000, LONG = FALSE)'
    } else if (method=='RE') {
      call <- 
        'Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2), OPT.TYPE = "RE", DATA = DATA, OPTNAME = "DE_12_RE", NUMITER = 200000, STATWINDOW = 50, DUMP.FREQ = 100000, LONG = FALSE)'
    }
    
  } else if (example==5) {
    text <- '
library(Optimus)
out.dir <- getwd()
region.len = 4e+04
mingap.Mb = 2e+06
out.name = "IJ.NEW.OPTI" 

seed = 840
opt.type = "SA" # "SA" | "RE"
# If opt.type = "RE" 
accratio = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2)
# If opt.type = "RE", replica should be equal to number of acceptance ratios 
# in accratio
replica = 4 # 4 | # 12
numiter = 2e+05 
cycles = 2
dump.freq = 1e+05
# For opt.type="SA", no need to set statwindow, default will be used (statwindow=70)
# For opt.type="RE", statwindow=50
statwindow = 50
calcErrorOnly = FALSE

#-------------------------------------------------------------------------------
# Function to calculate error rate of shuffled set of contacts in comparison to
# original set and based on the following criteria for a valid contact:
# 1. i should always be greater than j.  
# 2. i and j should be separated by a distance greater than 2 Mb. 
# 3. The set should only contain unique pairs/contacts.
# 4. Contact in new set should not be in the original set. 
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
# Define interfacing functions for Optimus
#-------------------------------------------------------------------------------
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

out.name <- paste0(out.name, ".", opt.type)

# Required gap between contacting bins; in terms of bins
mingap <- mingap.Mb/region.len

# Load IJ.ORIG
data("ij_orig")

if(calcErrorOnly==FALSE){
  
  DATA <- list(IJ.ORIG=IJ.ORIG, numContacts=nrow(IJ.ORIG), gaplimit=mingap)
  
  K <- IJ.ORIG
  
  set.seed(seed)
  # First random shuffle of js; is are kept in the same order
  K[,"j"] <- sample(x=K[,"j"], size=nrow(K), replace=FALSE)
  
  # Uses default SEED = 840
  
  if(opt.type=="SA"){
    Optimus(NCPU = replica, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = opt.type,
            OPTNAME = out.name, DATA = DATA, NUMITER = numiter, CYCLES = cycles, DUMP.FREQ = dump.freq,
            LONG = TRUE, SEED = seed)
  } else if(opt.type=="RE"){
    # Uses default SEED = 840
    Optimus(NCPU = replica, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = accratio,
            OPT.TYPE = opt.type, DATA = DATA, OPTNAME = out.name, NUMITER = numiter, STATWINDOW = statwindow,
            DUMP.FREQ = dump.freq, LONG = TRUE, SEED = seed)
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
'
    if (method=='SA') {
      call <- 
        ' Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA",OPTNAME = "IJ.NEW.OPTI.SA", DATA = DATA, NUMITER = 2e+05, CYCLES = 2, DUMP.FREQ = 1e+0 LONG = FALSE)'
    } else if (method=='RE') {
      call <- 
        'Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = ACCRATIO, OPT.TYPE = "RE", DATA = DATA, OPTNAME = "IJ.NEW.OPTI.RE", NUMITER = 2e+05, STATWINDOW = 50, DUMP.FREQ = 1e+05, LONG = FALSE)'
    }
  }
  fileConn<-file(file_name)
  writeLines(c(paste0('setwd("',dir,'/example_',example,'")'),text,call), fileConn)
  close(fileConn)
  if (run) {
    source(file_name)
  }
}