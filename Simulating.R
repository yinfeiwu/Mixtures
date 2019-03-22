# "snpsampgen.R" contains a function that takes a new set of 100,000 SNP and uses those to simulate a sample population of size k.

# "SLSQPmixturesR.R" contains a function that feeds this sample population matrix into the SLSQP algorithm in R.

# "HA_script.py" contains a function that feeds this sample population matrix into the SLSQP algorithm in Python.

# At the end of every iteration of this loop will be printed the estimated pi-values, number of iterations, and time for both versions of SLSQP.

# The source scripts needed to execute this code can be found at https://github.com/GregoryMatesi/SimulationDesign

# Thanks to Ian for "snpsampgen.R" and "SLSQPmixturesR.R." And to Jordan for "HA_script.py."

total.c <- read.table("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/5ancfinalsimdat.txt")
#install.packages("nloptr")    # Needs to be installed on the server. It contains our SLSQP R function.
#install.packages("reticulate")
library("nloptr")
library("reticulate")         # Reticulate is already installed on the server. It links Python with R.

source("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/snpsampgen.R")
source("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/SLSQPmixturesR.R")
source_python("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/Greg/HA_script.py")    # source_Python() from the "reticulate" package.
# push to Audrey

START <- Sys.time()
teamMember <- "Greg"              # Put your name here. Please make sure it matches the name of your folder in math-compute.
# setwd(pasteO("/nfs/storage/math/gr0ss-s2/projects/mixtures/team_members/current_team" , teamMember))
# Write the file to your personal folder in the math-compute server.

numberSims <- 3                        # Start with 10. then try 100 or 1000.
k <-5                                  # Keep at 5. This number is input for simulating the sampling population.
numAnc <- 2                            # Set by the user. This number is only for recording in the spreadsheet.
interval <- "(0.01_0.05)"              # Change your interval
ancestry1 <- "eur_MAF"                 # Name the ancestries you are working with
ancestry2 <- "afr_MAF"                 #

guess1 <- 1/k                          # Generally leave this along for now.           

output <- c()                          # Initialize an empty output matrix.



for (i in 1:numberSims){
    
    seed <- Sys.time() # Save the seed used in the uniform draw and for sampling population.
    
    set.seed(seed)
    eurfrac <- runif(1, 0.01, 0.05)    # Uniform random from user chosen interval.
    afrfrac <- 1 - eurfrac             
    easfrac <- 0
    sasfrac <- 0
    namfrac <- 0
    
    
    set.seed(seed)
    A <- snpsampgen(k, "eur_MAF" , "afr_MAF", "eas_MAF", "sas_MAF", "nam_MAF", eurfrac, afrfrac, easfrac, sasfrac, namfrac)    
    # Simulating a sample pop with new SNPs.
    
    # Calling the SLSQPmixtures function from SLSQPmixturesR.R
    Python_ <- print(SLSQPmixtures(A, k))                              # stores pi-values, iters, time.
    
    # Calling the HA function from HA_script.py                      
    af <- cbind(A$AF)                                                  # Total allele frequency vector
    A <- cbind(A$eur_MAF, A$afr_MAF, A$eas_MAF, A$sas_MAF, A$nam_MAF)  # CHANGEME: cbind(A$afr_MAF, A$CEU_MAF) etc
    
    guess <- rbind(guess1, guess1, guess1, guess1, guess1)             # kx1 vector of starting guesses.
    HA_ <- print(HA(A, af, guess))                                     # Stores pi-values, iters, time.
    
    output <- rbind(output, c(HA_[1], HA_[2], HA_[3], HA_[4], HA_[5], HA_[6], HA_[7], Python_[1], Python_[2], Python_[3], Python_[4], Python_[5], Python_[6], Python_[7], seed, teamMember, eurfrac, afrfrac, easfrac, sasfrac, namfrac))
    print(i)
    print(Sys.time() - START)
}
colnames(output) <- c("R.eur", "R.afr", "R.eas", "R.sas", "R.nam",  "R.iters", "R.time","Python.eur", "Python.afr", "Python.eas", "Python.sas", "Python.nam", "Python.iters", "Python.time", "seed", "Team.memeber", "True.eur", "True.afr", "True.eas", "True.sas", "True.nam")

# write output to personal file
# past0("/nfs/storage/math/gr0ss-s2/projects/mixtures/team_members/current_team/", teamMember, n, interval, ancestry1, "_", ancestry2)
# past0("/home/jovyan, teamMember, n, interval, ancestry1, "_", ancestry2)
#write.csv(output, file = paste0("/home/jovyan/", teamMember, "/", numAnc, "_", interval, "_", ancestry1, "_", ancestry2))
write.csv(output, file = paste0("/nfs/storage/math/gross-s2/projects/mixtures/team_members/current_team/", teamMember, "/", numAnc, "_", interval, "_", ancestry1, "_", ancestry2, ".csv"))


END <- Sys.time()
END - START
