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

Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", OPTNAME = "DE_4_SA", DATA = DATA, NUMITER = 200000, CYCLES = 2, DUMP.FREQ = 100000, LONG = FALSE)

Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2), OPT.TYPE = "RE", DATA = DATA, OPTNAME = "DE_12_RE", NUMITER = 200000, STATWINDOW = 50, DUMP.FREQ = 100000, LONG = FALSE)