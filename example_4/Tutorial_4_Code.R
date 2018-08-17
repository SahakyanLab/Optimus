# An object that is created by setting the state names and the initial values.
# The values are to be given in any unit. The sum will stay constant across the
# dynamics of the system. For instance, in case those are the fractions in %,
# the sum should be 100% and will stay the same in the solution as well.
state  <- c(cA=100, cB=100, cC=100, cAP=0, cBP=0, cCP=0)

# An object holding the target observables for the optimisation.
target <- c(cA=40, cB=20, cC=70, cAP=60, cBP=80, cCP=30)

# An object holding the model definition. It should contain equations that use
# the objects with the names specified within state and coef above, and should
# have equations that assign the outcomes to new objects that have the same
# order and names as specified in state, but with "d" at the beginning.
#---------------------------------
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

# An object holding all the coefficients used in the model. The setting of the
# names and the initialisation is done simultaneusly by editing this.

################################################################################
K <- c(k1=1.0, k2=1.0, k3=1.0, k4=1.0)
################################################################################

################################################################################
library(deSolve)
m <- function(K, DATA){
  state <- DATA$state
  model <- DATA$model
  
  span = 10.0  # byr

  times <- c(0, span)
  out   <- ode(y=state, times=times, func=model, parms=K)[2,2:(length(state)+1)]
  return(out)
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
  move.step <- 0.0005
  K.new[K.ind.toalter] <- K.new[K.ind.toalter] + sample(size=1, x=c(-move.step, move.step))
  
  ## Setting the negative coefficients to 0
  neg.ind <- which(K.new < 0)
  if(length(neg.ind)>0){ K.new[neg.ind] <- 0 }
  
  return(K.new)
}
################################################################################

Optimus(NCPU = 1, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", OPTNAME = "DE_4_SA", DATA = DATA, NUMITER = 100000, CYCLES = 2, DUMP.FREQ = 50000, LONG = FALSE)
