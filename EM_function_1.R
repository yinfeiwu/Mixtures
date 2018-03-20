####################################################################################
## code to estimate parameters of mixture distribution from frequency level data  ##
## 9/1/2017														                  ##
####################################################################################

##qsub  -q all.q submit_R.sh

########  THIS IS THE CODE THAT WORKS FOR AFRICAN AND EUROPEAN ONLY -- i.e. EM1 #########

setwd("/home/projects/mixtures")

#######  THIS DOES NOT WORK FOR AFRICAN AND EUROPEAN AND ALL SAMPLES -- i.e. EM2  #####

person="Kendra" ## make a directory with your name - make sure you use your name when running the code so you are not overwriting others work

########  sim parameters  ##########
Ntot=10000  ##total number of simulated samples
N_AFR=5000
MAF_thresh=0.05
pop_names<-c("AFR", "NFE") ##c("AFR","AMR", "EAS","FIN","OTH","SAS", "NFE") ##
k=length(pop_names)
pops<-ifelse(k==2, "AFR_NFE", "ALL") ##"ALL"  ##AFR_NFE"  ##
pi_start<-rep(1/k, k)  ##c(.5, .5) ##c(.9, .1)
threshold=0.0000001

#############  ExAC data  ###########
exac<-read.table("ExAC/all.common_added_cols_v2.txt", sep="\t", header=T, as.is=T, fill=T)


############# simulating NFE and African mixture ##########
x.obs.NFE<- as.data.frame(exac[,c("AC_hom_NFE", "AC_het_NFE", "AC_homref_NFE")])
x.obs.AFR<- as.data.frame(exac[,c("AC_hom_AFR", "AC_het_AFR", "AC_homref_AFR")])
x.sim.NFE<-t(apply(x.obs.NFE, 1, function(x){rmultinom(1, (Ntot-N_AFR), prob=(x))}))
x.sim.AFR<-t(apply(x.obs.AFR, 1, function(x){rmultinom(1, N_AFR, prob=(x))}))

x.touse.a<-x.sim.NFE+x.sim.AFR
x.touse<-as.data.frame(x.touse.a)
x.touse$AF<-(2*x.touse[,1]+x.touse[,2])/(2*Ntot)
x<-x.touse[x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh),1:3]

###########  EM1 FUNCTION  #############
N=nrow(x) # number of markers

if(pops=="AFR_NFE"){p<-exac[which(x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh)),c("AC_hom_AFR", "AC_het_AFR", "AC_homref_AFR","AC_hom_NFE", "AC_het_NFE", "AC_homref_NFE")]}
if(pops=="ALL"){p<-exac[which(x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh)),c("AC_hom_AFR", "AC_het_AFR", "AC_homref_AFR","AC_hom_AMR", "AC_het_AMR", "AC_homref_AMR","AC_hom_EAS", "AC_het_EAS", "AC_homref_EAS", "AC_hom_FIN", "AC_het_FIN", "AC_homref_FIN", "AC_hom_OTH", "AC_het_OTH", "AC_homref_OTH", "AC_hom_SAS", "AC_het_SAS", "AC_homref_SAS", "AC_hom_NFE", "AC_het_NFE", "AC_homref_NFE")]}

##p<-apply(p.a, 2, function(x){ifelse(is.na(x),0,x)})  ##assuming missing in ExAC is AF=0


pi_out<-pi_out_median<-pi_new<-pi_start
iter=0
thresh_check=threshold+1

while( thresh_check > threshold){
    pi=pi_new
    gamma_tmp.a<-matrix(NA, nrow=nrow(x), ncol=k)
    for(p_k in 1:k){
        ##estimating the posterior probs from the data per k ancestry group per n markers
        x.p<-cbind(x, p[,c(((p_k-1)*3+1):((p_k-1)*3+3))])
        gamma_tmp.a[,p_k]<-apply(x.p,1, function(x.f){
            p_k_geno<-as.numeric(x.f[4:6])
            pi[p_k]*dmultinom(as.numeric(x.f[1:3]), prob = p_k_geno, log=T)
        })
    }
    
    gamma_tmp<-(apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})-gamma_tmp.a)/apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})
   
    pi_new<-apply(gamma_tmp,2,function(x){mean(x, na.rm=T)})
    pi_median<-apply(gamma_tmp,2,function(x){median(x, na.rm=T)})
    pi_90<-apply(gamma_tmp,2,function(x){quantile(x, probs=0.90,na.rm=T)})
        
    pi_out<-rbind(pi_out, pi_new)
    pi_out_median<-rbind(pi_out_median, pi_median)
    iter=iter+1
    thresh_check<-sum(abs(pi-pi_new))
    if(round(iter/10)==(iter/10)){write.table(gamma_tmp, paste("gamma_values/pops", pops, "AFR", N_AFR, "_NFE", (N-N_AFR), "_start", paste(pi_start, collapse="_"), "_iter", iter, ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")}
    print(pi)
    print(nrow(x))
}


pi_out2<-cbind(pi_out, pi_out_median, iter, N)

write.table(pi_out2, paste(person, "/pops", pops, "_", N_AFR, "_", (Ntot-N_AFR), "_start", paste(pi_start, collapse="_"), ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")




