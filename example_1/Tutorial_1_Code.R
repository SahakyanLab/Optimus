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

Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", DATA = DATA, OPTNAME = "poly_4_SA", LONG = FALSE)

Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2), OPT.TYPE = "RE", DATA = DATA, OPTNAME = "poly_12_RE", LONG = FALSE)
