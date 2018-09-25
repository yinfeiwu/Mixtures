####################################################################################
## code to estimate parameters of mixture distribution from frequency level data  ##
## EM1 - precisely estimating mixture proportions for known ancestral groups      ##
## March 2018                                                                      ##
####################################################################################

##qsub  -q all.q submit_R.sh

setwd("/home/projects/mixtures")


#############  read in data  ###########
temp<-c()
for(chr in 15:22){
    load(paste("1000Genomes/global_populations/Global_ancestries_BROADfiltering.chr", chr,".RData", sep="")) ##global
    temp<-rbind(temp, global[!duplicated(global$pos),])
    print(chr)
}


temp2<-temp[sample(1:nrow(temp))[1:200000],]
temp<-temp2[temp2$AFR_AF>0.05 & temp2$EUR_AF>0.05 & temp2$AFR_AF<0.95 & temp2$EUR_AF<0.95, ]

##########  parameters  ###########
N_pop1<-2000
Ntot<-10000
pop_names<-c("EUR", "AFR")
MAF_thresh=0.05

############# simulating European and African mixture ##########
x.sim.pop2<-t(sapply(temp[,paste(pop_names[2], "_AF", sep="")], function(x){x2<-as.numeric(x); rmultinom(1, (Ntot-N_pop1), prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))

x.sim.pop1<-t(sapply(temp[,paste(pop_names[1], "_AF", sep="")], function(x){x2<-as.numeric(x); rmultinom(1, (N_pop1), prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))

x.touse.a<-x.sim.pop2+x.sim.pop1
x.touse<-data.frame(temp[,c(1:4)],x.touse.a)
x.touse$AF<-(2*x.touse[,5]+x.touse[,6])/(2*Ntot)
x<-x.touse[x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh),c(1,2,5:7)]

#########  EM data  ###########
dat<-temp[,c("chr", "pos", paste(pop_names[1], "_AF", sep=""), paste(pop_names[2], "_AF", sep=""))]

#########  EM parameters  ###########
threshold = 0.001
N <- nrow(x)
k=length(pop_names)
pi_init<-rep(1/k, k)
pi_out <- pi_out_median <- pi_out_90 <- pi_new <- pi_median <- pi_90 <- pi <- pi_init
names(pi_median) <- names(pi_new) <- names(pi_90) <- names(pi) <- pop_names
thresh_check <- threshold + 1

iter = 0
gamma_tmp.a.out<-gamma_tmp_new<-c()

###  checking that reference data and input data have the same number of observations  ###
x$chr_pos<-paste(x$chr, x$pos, sep="_")
dat$chr_pos<-paste(dat$chr, dat$pos, sep="_")

x2<-x[x$chr_pos %in% dat$chr_pos,]
dat2<-dat[dat$chr_pos %in% x$chr_pos,]

x<-x2[,c(1:5)]
dat<-dat2[,c(1:(2+k))]

##  check for if data is in genotypes  ##
if(ncol(dat)<(k*3+2)){dat2<-data.frame(dat[,1:2])
    for(i in 1:k){tmp<-as.numeric(dat[,(2+i)])
        dat2<-data.frame(dat2, tmp**2,tmp*(1-tmp)*2, (1-tmp)**2)}
    dat=dat2}

p <- dat # p <- dat(which(x$AF > MAF_thresh & x$AF < (1-MAF_thresh)))


###########  EM alogirtm  ###########
while( thresh_check > threshold){
    pi=pi_new
    gamma_tmp.a<-matrix(NA, nrow=nrow(x), ncol=k)
    for(p_k in 1:k){
        ##estimating the posterior probs from the data per k ancestry group per n markers
        x.p<-cbind(x, p[,c(((p_k-1)*3+3):((p_k-1)*3+5))])
        gamma_tmp.a[,p_k]<-apply(x.p,1, function(x.f){
            
            p_k_geno<-as.numeric(x.f[6:8])
            pi[p_k]*dmultinom(as.numeric(x.f[3:5]), prob = p_k_geno, log=T)
                    })
    }
    

    ## note, since we take the log above we need to subtract the probability from 1, which is what we do below ##
    gamma_tmp<- (apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})-gamma_tmp.a)/apply(gamma_tmp.a,1,function(x){sum(x, na.rm=T)})
    
    
       pi_new<-apply(gamma_tmp,2,function(x){mean(x, na.rm=T)})
    pi_median<-apply(gamma_tmp,2,function(x){median(x, na.rm=T)})
    
    pi_out<-rbind(pi_out, pi_new)
    pi_out_median<-rbind(pi_out_median, pi_median)
    iter=iter+1
    if(iter>1){thresh_check<-sum(abs(gamma_tmp_new-gamma_tmp), na.rm=T)}
    gamma_tmp_new<-gamma_tmp
    thresh_check<-sum(abs(pi-pi_new))
    ##if(round(iter/10)==(iter/10)){write.table(gamma_tmp, paste(person, "/gamma_values/pops", pop_names[1], N_pop1, "_", pop_names[2], (N-N_pop1), "_start", paste(pi_start, collapse="_"), "_iter", iter, ".txt",sep=""), row.names=F, col.names=F, quote=F, sep="\t")}
    print(pi)
    print(nrow(x))
    ##print(table(is.na(gamma_tmp[,1])))
}

pi ##mean values
pi_median  ##median values


