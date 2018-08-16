# An object holding all the coefficients used in the model. The setting of the
# names and the initialisation is done simultaneusly by editing this.
# coef is REQUIRED by mcara()!
coef <- c(i=1.0, j=1.0, k=1.0,
          l=1.0, m=1.0, n=1.0)
          

# An object holding the target observables for the optimisation.
# target is REQUIRED by mcara()!
target <- c(Ca=29.7548, Cg=20.23669, Ct=29.793, Cc=20.21551)



# An object that is created by setting the state names and the initial values.
# The values are to be given in any unit. The sum will stay constant across the
# dynamics of the system. For instance, in case those are the fractions in %,
# the sum should be 100% and will stay the same in the solution as well.
state  <- c(Ca=100/4, Cg=100/4, Ct=100/4, Cc=100/4)

# An object holding the model definition. It should contain equations that use
# the objects with the names specified within state and coef above, and should
# have equations that assign the outcomes to new objects that have the same
# order and names as specified in state, but with "d" at the beginning.
#---------------------------------
model <- function(t, state, coef){

  with( as.list(c(state, coef)), {
    # rate of change
       dCa <- i*Cc+l*Ct+m*Cg-(j+l+n)*Ca
       dCg <- n*Ca+k*Cc+j*Ct-(m+i+k)*Cg
       dCt <- l*Ca+i*Cg+m*Cc-(l+n+j)*Ct
       dCc <- j*Ca+n*Ct+k*Cg-(i+m+k)*Cc
    # return the rate of change
    list(c(dCa, dCg, dCt, dCc))
  }) # end with(as.list ...

}
#---------------------------------

# Helper function.
#---------------------------------
GetModelOut <- function(state=state, model=model, coef=coef){
  span = 10.0  # byr
  library("deSolve")

  times <- c(0, span)
  out   <- ode(y=state, times=times, func=model, parms=coef)[2,2:(length(state)+1)]
  return(out)
}
#---------------------------------

# Helper function to assess the agreement quality.
#---------------------------------
GetPseudoE <- function(sol=sol, target=target, type="RMSE"){
  RESULT <- NULL
  if(type=="RMSE"){
    RESULT$Q <- sqrt(mean((sol-target)^2)) # measure of the fit quality
    RESULT$E <- RESULT$Q # the pseudoenergy derived from the above measure
  }
  return(RESULT)
}
#---------------------------------


# The modelling core function that is called in the main mcara() optimiser.
# Takes - coef and other arguments specified in this file (model.R), outputs -
# an object with $sol, $sol.perf and $coef.new components required by mcara.R.
#---------------------------------
MDLCore <- function(coef=coef, move.step=0.005, state=state,
                    model=model, target=target){                               
  coef.new <- coef
  coef.ind.toalter <- sample(size=1, x=1:length(coef.new))
  # Creating a potentially new set of coefficients where one entry is altered
  # by either +move.step or -move.step:
  coef.new[coef.ind.toalter] <- coef.new[coef.ind.toalter] +
                                sample(size=1, x=c(-move.step, move.step))
  sol      <- GetModelOut(state=state, model=model, coef=coef.new)
  sol.perf <- GetPseudoE(sol=sol, target=target, type="RMSE")
  
  RESULT          <- NULL
  RESULT$sol      <- sol
  RESULT$coef.new <- coef.new
  RESULT$sol.perf <- sol.perf
  return(RESULT)

}
#---------------------------------


################################################################################
