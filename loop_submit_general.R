##########################
## loop submit command  ##
## Audrey Hendricks     ##
## April 3 2018         ##
## Updated by Megan     ##
## April 10 2018
##########################

##########################
##########################
## Make sure you are in your loop directory!!!
##########################
##########################

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
  
  write(c(paste("Ntot = ", Ntot, sep=""), paste("N_pop1 = ", N_pop1, sep=""), paste("MAF_thresh = ", MAF_thresh, sep=""), paste("pop_names = c('", paste(t(pop_names_loop[i,]), collapse="', '"), "')", sep=""), paste("threshold = ", threshold, sep="")), file=paste(tmp, ".R", sep=""))
  write("source('/home/projects/mixtures/EM_function_1_forlooping.R')", file=paste(tmp, ".R", sep=""), append=T)
  
  write(c("#$ -cwd",paste("#$ -o ", getwd(), '/', tmp, ".log", sep=""),paste("#$ -o ",getwd(), '/', tmp, ".err", sep=""),"#$ -S /bin/bash", "", paste("Rscript ", getwd(), '/', tmp, ".R", sep="")),file=paste(tmp, ".sh", sep=""))
  
  tosubmit<-c(tosubmit, paste(getwd(), '/', tmp, ".sh", sep=""))
}

write(tosubmit, 'tosubmit.txt')


