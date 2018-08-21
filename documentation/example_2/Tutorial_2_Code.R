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
  E <- AIC(O)/1000 # Akaike's information criterion.

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

Optimus(NCPU = 4, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, OPT.TYPE = "SA", DATA = DATA, OPTNAME = "term_4_SA", NUMITER = 200000, CYCLES = 2, DUMP.FREQ = 100000, LONG = FALSE)

Optimus(NCPU = 12, K.INITIAL = K, rDEF = r, mDEF = m, uDEF = u, ACCRATIO = c(90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10, 2), OPT.TYPE = "RE", DATA = DATA, OPTNAME = "term_12_RE", NUMITER = 200000, STATWINDOW = 50, DUMP.FREQ = 100000, LONG = FALSE)
