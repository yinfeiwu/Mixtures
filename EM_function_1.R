####################################################################################
## code to estimate parameters of mixture distribution from frequency level data  ##
## EM1 - precisely estimating mixture proportions for known ancestral groups      ##
## March 2018														              ##
####################################################################################

##qsub  -q all.q submit_R.sh

setwd("/home/projects/mixtures")

#######  THIS DOES NOT WORK FOR unknown reference sampels -- i.e. EM2  #####

person="Kendra" ## make a directory with your name - make sure you use your name when running the code so you are not overwriting others work


#############  ExAC data  ###########
exac<-read.table("ExAC/all.common_added_cols_v2.txt", sep="\t", header=T, as.is=T, fill=T)


########  sim parameters  ##########
k=length(pop_names)


############# simulating NFE and African mixture ##########
x.obs.pop2<- as.data.frame(exac[,paste(c("AC_hom_", "AC_het_", "AC_homref_"), pop_names[2], sep="")])
x.obs.pop1<- as.data.frame(exac[,paste(c("AC_hom_", "AC_het_", "AC_homref_"), pop_names[1], sep="")])
x.sim.pop2<-t(apply(x.obs.pop2, 1, function(x){rmultinom(1, (Ntot-N_pop1), prob=(x))}))
x.sim.pop1<-t(apply(x.obs.pop1, 1, function(x){rmultinom(1, N_pop1, prob=(x))}))

x.touse.a<-x.sim.pop2+x.sim.pop1
x.touse<-as.data.frame(x.touse.a)
x.touse$AF<-(2*x.touse[,1]+x.touse[,2])/(2*Ntot)
x<-x.touse[x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh),1:3]

###########  EM1 FUNCTION  #############
N=nrow(x) # number of markers

p<-exac[which(x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh)),paste(c("AC_hom_", "AC_het_", "AC_homref_"), rep(pop_names, each=3), sep="")]


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
    if(round(iter/10)==(iter/10)){write.table(gamma_tmp, paste(person, "/gamma_values/pops", pop_names[1], N_pop1, "_", pop_names[2], (N-N_pop1), "_start", paste(pi_start, collapse="_"), "_iter", iter, ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")}
    print(pi)
    print(nrow(x))
}


pi_out2<-cbind(pi_out, pi_out_median, iter, N)

write.table(pi_out2, paste(person, "/pops", pop_names[1], "_", N_pop1, "_", pop_names[2], (Ntot-N_pop1), "_start", paste(pi_start, collapse="_"), ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")




