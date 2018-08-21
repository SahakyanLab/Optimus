################################################################################
K <- c(PHI=90, PSI=90)
################################################################################

################################################################################
m <- function(K,
              notconvergedE = -100.00,
              mopac.cmd="/home/alex/prog/mopac2016/MOPAC2016.exe",
              clean = TRUE){

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
