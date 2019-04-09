# "snpsampgen.R" contains a function that takes a new set of 100,000 SNP and uses those to simulate a sample population of size k.

# "SLSQPmixturesR.R" contains a function that feeds this sample population matrix into the SLSQP algorithm in R.

# "HA_script.py" contains a function that feeds this sample population matrix into the SLSQP algorithm in Python.

# At the end of every iteration of this loop will be printed the estimated pi-values, number of iterations, and time for both versions of SLSQP.

# The source scripts needed to execute this code can be found at https://github.com/GregoryMatesi/SimulationDesign

# Thanks to Ian for "snpsampgen.R" and "SLSQPmixturesR.R." And to Jordan for "HA_script.py."

total.c <- read.table("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/5ancfinalsimdat.txt")

library("nloptr")
library("reticulate")         # Reticulate is already installed on the server. It links Python with R.

source("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/snpsampgen.R")
source("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/SLSQPmixturesR.R")
source_python("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/HA_script.py")    # source_Python() from the "reticulate" package.


START <- Sys.time()
teamMember <- "Greg"              # Put your name here. Please make sure it matches the name of your folder in math-compute.
interval1 <- 0.015                # Change this to the interval you are working on
interval2 <- 0.025                # 

numberSims <- 3                        # Start with 3. then try 10. Then 1000!
k <- 5                                 # Keep at 5. This number is input for simulating the sampling population.
numAnc <- 2                            # Set by the user. This number is only for recording in the spreadsheet.
ancestry1 <- "afr_MAF"                 # Name the ancestries you are working with
ancestry2 <- "nam_MAF"                 #

guess1 <- 1/k                          # Generally leave this alone for now.           

output <- c()                          # Initialize an empty output matrix.



for (i in 1:numberSims){
    
    seed <- Sys.time() # Save the seed used in the uniform draw and for sampling population.
    
    set.seed(seed)
    afrfrac <- runif(1, interval1, interval2) # Uniform random from user chosen interval.
    easfrac <- 0    
    sasfrac <- 0
    eurfrac <- 0
    namfrac <- 1 - afrfrac
    
    
    set.seed(seed)
    A <- snpsampgen(k, "afr_MAF", "eas_MAF", "sas_MAF", "eur_MAF" , "nam_MAF", afrfrac, easfrac, sasfrac, eurfrac, namfrac)    
    # Simulating a sample population with a new set of 100,00 SNPs.
    
    # Calling the SLSQPmixtures function from SLSQPmixturesR.R
    Python_ <- print(SLSQPmixtures(A, k))                              # Returns pi-values, iterations, time.
    
    # Calling the HA function from HA_script.py                      
    af <- cbind(A$AF)                                                  # Total allele frequency vector.
    A <- cbind(A$afr_MAF, A$eas_MAF, A$sas_MAF, A$eur_MAF, A$nam_MAF)  # All global populations used.
    
    guess <- rbind(guess1, guess1, guess1, guess1, guess1)             # kx1 vector of starting guesses.
    HA_ <- print(HA(A, af, guess))                                     # Returns pi-values, iterations, time.
    
    output <- rbind(output, c(HA_[1], HA_[2], HA_[3], HA_[4], HA_[5], HA_[6], HA_[7], Python_[1], Python_[2], Python_[3], Python_[4], Python_[5], Python_[6], Python_[7], seed, teamMember, afrfrac, easfrac, sasfrac, eurfrac, namfrac, eurfrac - as.numeric(HA_[1]), afrfrac - as.numeric(HA_[2]), easfrac - as.numeric(HA_[3]), sasfrac - as.numeric(HA_[4]), namfrac - as.numeric(HA_[5]), eurfrac - as.numeric(Python_[1]), afrfrac - as.numeric(Python_[2]), easfrac - as.numeric(Python_[3]), sasfrac - as.numeric(Python_[4]), namfrac - as.numeric(Python_[5]), interval1, interval2))
    print(i)
    print(Sys.time() - START)
}
colnames(output) <- c("R.afr", "R.eas", "R.sas","R.eur", "R.nam",  "R.iters", "R.time", "Python.afr", "Python.eas", "Python.sas", "Python.eur", "Python.nam", "Python.iters", "Python.time", "seed", "Team.memeber", "True.afr", "True.eas", "True.sas", "True.eur", "True.nam", "accuracy.HA.afr", "accuracy.HA.eas", "accuracy.HA.sas", "accuracy.HA.eur", "accuracy.HA.nam", "accuracy.Python.afr", "accuracy.Python.eas", "accuracy.Python.sas", "accuracy.Python.eur", "accuracy.Python.nam", "interval.1", "interval.2")

write.csv(output, file = paste0("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/", teamMember, "/", ancestry1, "_", interval1, "-", interval2, "_", ancestry2,  ".csv"))


END <- Sys.time()
END - START
