####################################################################################
## code to estimate parameters of mixture distribution from frequency level data  ##
## 9/1/2017														                  ##
####################################################################################

setwd("/Users/audreyhendricks/Desktop/CU-Denver/research/Active/mixture_Exac/")
#############  chr 20  ###########
chr20.common<-read.table("chr20.common.txt", sep="\t", header=T, as.is=T)

chr20.common$AC_het=sapply(chr20.common$INFO, function(x){unlist(strsplit(unlist(strsplit(x, "AC_Het="))[2], ";"))[1]})
chr20.common$AC_hom=sapply(chr20.common$INFO, function(x){unlist(strsplit(unlist(strsplit(x, "AC_Hom="))[2], ";"))[1]})
chr20.common$AC_homref<-chr20.common$an/2-as.numeric(chr20.common$AC_het)-as.numeric(chr20.common$AC_hom)


## five possible groups ##
x.obs<- as.data.frame(chr20.common[,c("AC_hom", "AC_het", "AC_homref")])  ##observed markers n x

## simulating a sample from the overall allele counts ##
## with noise ##
x.sim<-t(apply(x.obs, 1, function(x){rmultinom(1, 1000, prob=(t(as.numeric(x))*t(runif(3, min=0.8, max=1.2))))}))
## without noise ##
x.sim<-t(apply(x.obs, 1, function(x){rmultinom(1, 1000, prob=(x))}))
x<-x.sim

N=nrow(x) # number of markers
p.a<- chr20.common[,c("AF_AFR", "AF_AMR", "AF_EAS", "AF_FIN", "AF_OTH", "AF_SAS", "AF_NFE")] ## backbone probabilities -- from 1000Genomes or other reference data
p<-apply(p.a, 2, function(x){ifelse(is.na(x),0,x)})  ##assuming missing in ExAC is AF=0
k=ncol(p)  ## number of groups
z<-rep(1/k,k)  #z - group assignment##maf of each marker in each ancestral group; n by k matrix
##z<-c(0,0,0,0,1,0,0) ## can try different starting vectors of Z; seems to be stable
pi<-pi_new<- rep(1/k,k) ##vector of proportion of each ancestral group -- WHAT WE WANT TO ESTIMATE -- initialize as 1/k vector of k length


threshold=0.01 #threshold of precision to stop iterations
thresh_check=1 #initializing the threshold check value
iter.count=0 #initializing counting the number of iteractions

while( thresh_check > threshold){  ## 4 - check convergence
pi=pi_new	

gamma_tmp.a<-matrix(NA, nrow=N, ncol=k)
for(p_k in 1:k){

##estimating the posterior probs from the data per k ancestry group per n markers
	x.p<-cbind(x, p[,p_k])
gamma_tmp.a[,p_k]<-apply(x.p,1, function(x.f){
	p_k_l<-as.numeric(x.f[4])
p_k_geno<-matrix(c(p_k_l**2, 2*p_k_l*(1-p_k_l), (1-p_k_l)**2), ncol=3)  ##n by 3
dmultinom(as.numeric(unlist(x.f[1:3])), prob = unlist(p_k_geno))})
}

## calculating the probabilities per k ancestry group per n markers weighting by total prob values across all k ancestry groups -- estimate of the proportion of the marker that is from each ancestry group
gamma_tmp<-gamma_tmp.a/apply(gamma_tmp.a,1,sum)

## summing probabilities across all markers to arrive at mixture proportions and dividing by sum to ensure proportions sum to 1
N_k = apply(gamma_tmp, 2, function(x){sum(x, na.rm=T)})  ## vector of N for each group K
pi_new.a=N_k/N  ## vector of k length
pi_new<-pi_new.a/sum(pi_new.a)

## checking threshold and increasing iter count
## right now I am treating all groups the same, but we could weight distance away from previous proportion by expected sample size, etc.
thresh_check<-sqrt(sum((pi - pi_new)**2))
iter.count=iter.count+1
}## end while loop

