#####################################################################################
##  code to estimate parameters of mixture distribution from frequency level data  ##
##  EM2 - estimating hidden ancestries                                             ##
##  March 2018                                                                     ##
####################################################################################

##qsub  -q all.q submit_R.sh

########  EM2 ALGORITHM #########

setwd("/home/projects/mixtures")

#######  THIS DOES NOT estimate precise proportions -- i.e. EM1  #####

person="Kendra" ## make a directory with your name - make sure you use your name when running the code so you are not overwriting others work

#############  ExAC data  ###########
exac<-read.table("ExAC/all.common_added_cols_v2.txt", sep="\t", header=T, as.is=T, fill=T)


########  sim parameters  ##########
N_pop=c(5000,5000)  ## the order here should match the order of the sim.pop_names below
Ntot=sum(N_pop)
sim.pop_names<-c("AFR", "NFE")
MAF_thresh=0.01

ref.pop_names<-c("AFR","AMR", "EAS","FIN","OTH","SAS", "NFE")
k=length(ref.pop_names)
pi_start<-rep(1/k, k)
threshold=0.01  ## may need to depend on number of reference populations



############# simulating NFE and African mixture ##########
x.touse.a<-matrix(0,nrow=nrow(exac), ncol=3)

for(i in 1:length(sim.pop_names)){
x.obs<- as.data.frame(exac[,paste(c("AC_hom_", "AC_het_", "AC_homref_"), sim.pop_names[i], sep="")])
x.sim<-t(apply(x.obs, 1, function(x){rmultinom(1, N_pop[i], prob=(x))}))

x.touse.a<-x.touse.a+x.sim
}

x.touse<-as.data.frame(x.touse.a)
x.touse$AF<-(2*x.touse[,1]+x.touse[,2])/(2*Ntot)
x<-x.touse[x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh),1:3]

###########  FUNCTION  #############
N=nrow(x) # number of markers

p<-exac[which(x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh)),c("AC_hom_AFR", "AC_het_AFR", "AC_homref_AFR","AC_hom_AMR", "AC_het_AMR", "AC_homref_AMR","AC_hom_EAS", "AC_het_EAS", "AC_homref_EAS", "AC_hom_FIN", "AC_het_FIN", "AC_homref_FIN", "AC_hom_OTH", "AC_het_OTH", "AC_homref_OTH", "AC_hom_SAS", "AC_het_SAS", "AC_homref_SAS", "AC_hom_NFE", "AC_het_NFE", "AC_homref_NFE")]



pi_out<-pi_out_median<-pi_new<-pi_start
iter=0
thresh_check=threshold+1
filtering.thresh=.000001

print(pi_new)
print(nrow(x))

while( thresh_check > threshold){
    pi=pi_new
    gamma_tmp.a<-matrix(NA, nrow=nrow(x), ncol=k)
    for(p_k in 1:k){
        ##estimating the posterior probs from the data per k ancestry group per n markers
        x.p<-cbind(x, p[,c(((p_k-1)*3+1):((p_k-1)*3+3))])
        gamma_tmp.a[,p_k]<-apply(x.p,1, function(x.f){
            p_k_geno<-as.numeric(x.f[4:6])
            pi[p_k]*dmultinom(as.numeric(x.f[1:3]), prob = p_k_geno, log=F)
               })
    }
    
    gamma_tmp.c<-gamma_tmp.a[which(apply(gamma_tmp.a, 1, sum)>(filtering.thresh*iter)),]
    x<-x[which(apply(gamma_tmp.a, 1, sum)>(filtering.thresh*iter)),]
    p<-p[which(apply(gamma_tmp.a, 1, sum)>(filtering.thresh*iter)),]
    
    
    gamma_tmp<-(gamma_tmp.c)/apply(gamma_tmp.c,1,function(x){sum(x, na.rm=T)})

    
    pi_new<-apply(gamma_tmp,2,function(x){mean(x, na.rm=T)})
    pi_median<-apply(gamma_tmp,2,function(x){median(x, na.rm=T)})
    
    pi_out<-rbind(pi_out, pi_new)
    pi_out_median<-rbind(pi_out_median, pi_median)
    iter=iter+1
    thresh_check<-sum(abs(pi-pi_new))
    ##if(round(iter/10)==(iter/10)){write.table(gamma_tmp, paste(person, "/gamma_values/EM2_simpops", paste(paste(sim.pop_names, N_pop, sep=""), collapse="_"), "_startpops", paste(paste(ref.pop_names,round(pi_start,digits=2), sep=""), collapse="_"), "_iter", iter, ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")}
    print(pi_new)
    print(nrow(x))
}


pi_out2<-cbind(pi_out, pi_out_median, iter, N)

write.table(pi_out2, paste(person, "/EM2_simpops", paste(paste(sim.pop_names, N_pop, sep=""), collapse="_"), "_startpops", paste(paste(ref.pop_names,round(pi_start,digits=2), sep=""), collapse="_"), ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")




