##########################
## loop submit command  ##
## Audrey Hendricks     ##
## April 3 2018         ##
##########################

#### example of how to create many files  to submit at once. Try to submit just a few files the first time to make sure you are submitting the jobs correctly.  ####

#### NOTE! the original R code parameters

## create folder called "loop" to store your looped output ##
## in your own folder ##
mkdir("loop")




###########  create R and sh scripts to submit  ############
## run in R ##

########  sim parameters loop  ##########
## parameters that you want to loop over ##
## always check your paste code ##

pop_names_loop<-data.frame(pop1="AFR", pop2=c("AMR", "EAS","FIN","OTH","SAS", "NFE"))


########  other sim parameters  ##########
Ntot=10000  ##total number of simulated samples
N_pop1=5000
MAF_thresh=0.05
##pop_names<-c("AFR", "NFE") ##c("AFR","AMR", "EAS","FIN","OTH","SAS", "NFE") ##
threshold=0.0001

tosubmit<-c()
for(i in 1:nrow(pop_names_loop)){
    ###FIX###
   tmp<-paste("EM1_loop_Ntot", Ntot, "_Npop1", N_pop1, "_MAFthresh", MAF_thresh, "_popnames", paste(t(pop_names_loop[i,]), collapse=""), "_threshold", threshold, sep="")
    
    write(c(paste("Ntot = ", Ntot, sep=""), paste("N_pop1 = ", N_pop1, sep=""), paste("MAF_thresh = ", MAF_thresh, sep=""), paste("pop_names = c('", paste(t(pop_names_loop[i,]), collapse="', '"), "')", sep=""), paste("threshold = ", threshold, sep="")), file=paste("/home/projects/mixtures/loop/", tmp, ".R", sep=""))
    write("source('/home/projects/mixtures/EM_function_1_forlooping.R')", file=paste("/home/projects/mixtures/loop/", tmp, ".R", sep=""), append=T)
    
    write(c("#$ -cwd",paste("#$ -o /home/projects/mixtures/loop/", tmp, ".log", sep=""),paste("#$ -o /home/projects/mixtures/loop/", tmp, ".err", sep=""),"#$ -S /bin/bash", "", paste("Rscript /home/projects/mixtures/loop/", tmp, ".R", sep="")),file=paste("/home/projects/mixtures/loop/", tmp, ".sh", sep=""))
    
    tosubmit<-c(tosubmit, paste("/home/projects/mixtures/loop/", tmp, ".sh", sep=""))
}

write(tosubmit, "/home/projects/mixtures/loop/tosubmit.txt")

########  loop to submit R scripts  #########
## run in BASH  ##
for i in $(cat tosubmit.txt)
do
qsub -q all.q  $i
echo $i
done

