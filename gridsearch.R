
###############################
##  1000Genomes gridsearch   ##
##  9/25/2018                ##
###############################




setwd("/home/projects/mixtures/1000Genomes/global_populations")

temp<-c()
for(chr in 10:22){
    load(paste("Global_ancestries_BROADfiltering.chr", chr,".RData", sep="")) ##global
    temp<-rbind(temp, global[!duplicated(global$pos),])
    print(chr)
}

temp2<-temp[sample(1:nrow(temp))[1:200000],]
temp<-temp2[temp2$AFR_AF>0.05 & temp2$EUR_AF>0.05 & temp2$AFR_AF<0.95 & temp2$EUR_AF<0.95, ]

############# simulating European and African mixture ##########
N_pop1<-200
Ntot<-10000
pop_names<-c("EUR", "AFR")
MAF_thresh=0.05

x.sim.pop2<-t(sapply(temp[,paste(pop_names[2], "_AF", sep="")], function(x){x2<-as.numeric(x); rmultinom(1, (Ntot-N_pop1), prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))

x.sim.pop1<-t(sapply(temp[,paste(pop_names[1], "_AF", sep="")], function(x){x2<-as.numeric(x); rmultinom(1, (N_pop1), prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))

x.touse.a<-x.sim.pop2+x.sim.pop1
x.touse<-data.frame(temp[,c(1:4)],x.touse.a)
x.touse$AF<-(2*x.touse[,5]+x.touse[,6])/(2*Ntot)
x<-x.touse[x.touse$AF>MAF_thresh & x.touse$AF<(1-MAF_thresh),c(1,2,5:7)]



###########  grid search  #########
x_af<-x[,1:2]
x_af$af<-apply(x[3:5], 1, function(x.tmp){(x.tmp[1]*2+x.tmp[2])/(sum(x.tmp)*2)})



###  for 3 ancestries  ###
p_af<-temp[,c("chr", "pos", "EUR_AF", "AFR_AF", "EAS_AF")] ##limit to only 3 ancestries
x.p<-data.matrix(merge(x_af, p_af, by=c("chr", "pos"))) ##cbind(x_af, p_af[,-c(1:2)])

output<-c()
pi.a<-seq(0,1, by=0.01)

for(j in pi.a){
    for(i in pi.a){
        pi.tmp<-c(i, j, (1-i-j))
        tmp2<-t(pi.tmp %*% t(x.p[,-c(1:3)]))
        least.sq<-sum((tmp2-x.p[,"af"])**2)
        deviance<-sum(abs(tmp2-x.p[,"af"]))
        ll = logLik(glm(x.p[,"af"]~tmp2))
        output = rbind(output,c(ll,least.sq, deviance,i,j,1-i-j))
    }}


output[output[,2]==min(output[,2]),]
output[output[,3]==min(output[,3]),]
output[output[,1]==max(output[,1]),]

########  2 ancestries ###########
x_af<-x[,1:2]
x_af$af<-apply(x[3:5], 1, function(x.tmp){(x.tmp[1]*2+x.tmp[2])/(sum(x.tmp)*2)})

p_af<-temp[,c("chr", "pos", "EUR_AF", "AFR_AF")] ##limit to only 3 ancestries
x.p<-data.matrix(merge(x_af, p_af, by=c("chr", "pos"))) ##cbind(x_af, p_af[,-c(1:2)])

output<-c()
pi.a<-seq(0,1, by=0.01)


for(i in pi.a){
    pi.tmp<-c(i, (1-i))
    tmp2<-t(pi.tmp %*% t(x.p[,-c(1:3)]))
    least.sq<-sum((tmp2-x.p[,"af"])**2)
    deviance<-sum(abs(tmp2-x.p[,"af"]))
    ll = logLik(glm(x.p[,"af"]~tmp2))
    output = rbind(output,c(ll,least.sq, deviance,i,1-i))
}

output[output[,2]==min(output[,2]),]
output[output[,3]==min(output[,3]),]
output[output[,1]==max(output[,1]),]





