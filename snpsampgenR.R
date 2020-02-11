#load main data
#MAKE SURE DIRECTORY IS SET CORRECTLY
#load(file="mixtures/Mixtures.git/total_strict_fine_maf01_atleastone.RData")

# TODO: generate and save random seed, save in final table

snpsampgen <- function(k, ancestry1, ancestry2, popfrac1, popfrac2){
    snpsamps <- total.c[sample(1:nrow(total.c))[1:100000],]


    # total number of populations to include and ratios for given pops
    # list as vectors/arrays, in order
    # eg c('afr_MAF','eur_MAF','eas_MAF') c(1/3,1/5,7/15) c(.2 .3 .2)

    N_total = 10000 # total number of samples (people)
    MAF_thresh <- 0.01 # keep at 0.01
    pop_number = k # number of total ancestries in pop, MAKE SURE pop_number >= length(pop_list)
    pop_list = c(ancestry1, ancestry2) # exact column names, complete random -> c()
    pop_fractions = c(popfrac1, popfrac2) # sum =< 1, final number disregarded if exact pop list, complete random -> c()


    # create table with all ancestry names and column locations
    # pulls any columns ending with '_MAF'

    fsdindices <- (grep("_MAF", colnames(total.c)))
    namevar <- colnames(total.c)
    fsdnames <- numeric(length(fsdindices))
    for (i in 1:length(fsdindices)){
      fsdnames[i] <- namevar[fsdindices[i]]
    }
    mafind <- data.frame(fsdindices, fsdnames)
    rm(fsdindices,namevar,fsdnames)

    # create table with correct population numbers
    # fill in missing values with randomized values

    # generate ancestry list
    tot_pop_list <- numeric(pop_number)
    input_pop <- which(mafind[,2] %in% pop_list)
    if(length(input_pop) > 0){ # adds user ancestry list
      for (i in 1:length(input_pop)){
        tot_pop_list[i] <- pop_list[i]
      }
    }
    if(length(input_pop) < pop_number){ # generates random ancestries
      rand_pops_add <- sample(mafind$fsdnames[!(mafind$fsdnames %in% pop_list)], (pop_number - length(input_pop)), replace = FALSE)
      tot_pop_list[(1+length(pop_list)):pop_number] <- as.character(rand_pops_add)
      rm(rand_pops_add)
    }

    # generate fraction list, sum to 1
    tot_frac_list <- numeric(pop_number)
    if (length(pop_fractions) > 0){ # adds user fraction list
      tot_frac_list[1:length(pop_fractions)] <- pop_fractions
    }
    if (sum(tot_frac_list) < 1 && !(length(pop_fractions) == pop_number)){ # generates random fractions if needed
      for (i in (1 + length(pop_fractions)):pop_number){
        if (i == pop_number){
          tot_frac_list[i] <- (1 - sum(tot_frac_list))
        }
        else {
          tot_frac_list[i] <- runif(1, min = .01, max = ((1 - sum(tot_frac_list)) - (.01 * (pop_number - i)))) # tweak to set min ancestry pop, MIN MUST EQUAL (TOT OPEN VALS LEFT * MIN)
        }
      }
    }

    # generates indices list for user/random ancestries
    tot_ind_list <- numeric(pop_number)
    for (i in 1:pop_number){
      tot_ind_list[i] <- mafind$fsdindices[mafind$fsdnames == tot_pop_list[i]]
    }

    # generates populations for user/random ancestries
    tot_popnum_list <- numeric(pop_number)
    for (i in 1:pop_number){
      if (i == pop_number){
        tot_popnum_list[i] <- (N_total - sum(tot_popnum_list))
      }
      else {
        tot_popnum_list[i] <- floor((N_total * tot_frac_list[i]))
      }
    }

    # combine into final table
    pop_frame <- data.frame(tot_ind_list,tot_pop_list,tot_frac_list,tot_popnum_list)
    rm(tot_frac_list,tot_ind_list,tot_pop_list,tot_popnum_list,i,input_pop)

    # ERROR WARNING
    # SOME PROB VALUES ARE NA VALUES
    # Causing error : Error in rmultinom(1, pop_frame$tot_popnum_list[i], prob = (c(x2^2, 2 *  : 
    # NA in probability 
    # TODO: filter possible NA values, potential solution resample SNPs, potential solution set NA to 0
    pop_matrix <- matrix(0, nrow = 100000, ncol = 3)
    for (i in 1:pop_number){
      popmatrixadd <- t(sapply(snpsamps[[pop_frame$tot_ind_list[i]]], function(x){x2<-as.numeric(x); rmultinom(1, pop_frame$tot_popnum_list[i], prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
      pop_matrix <- pop_matrix + popmatrixadd
    }
    rm(popmatrixadd)

    # Creating allele frequencies from simulated data
    # copy pasted from Audreys code

    master_frame_gen1 <- data.frame(snpsamps[,c(1:4)], pop_matrix)
    master_frame_gen1$AF <- (2 * master_frame_gen1[,5] + master_frame_gen1[,6]) / (2 * N_total)
    master_frame_gen2 <- master_frame_gen1[master_frame_gen1$AF > MAF_thresh & master_frame_gen1$AF < (1-MAF_thresh), c(1,2,5:7)]

    master_frame_AF <- master_frame_gen2[,1:2]
    master_frame_AF$AF <- apply(master_frame_gen2[3:5], 1, function(master_frame_gen2.tmp){(master_frame_gen2.tmp[1] * 2 + master_frame_gen2.tmp[2])/(sum(master_frame_gen2.tmp) * 2)})

    snpsamps_AF <- snpsamps[snpsamps$SNP %in% master_frame_AF$SNP,]

    samp_af_list <- numeric(2 + pop_number)
    samp_af_list[1:2] <- c(1,2)
    for (i in 3:(2 + pop_number)){
      samp_af_list[i] <- pop_frame$tot_ind_list[(i-2)]
    }
    snpsamps_AF <- snpsamps_AF[,samp_af_list]

    master_frame_final <- merge(snpsamps_AF, master_frame_AF, by = c('CHR', 'SNP'))
    rm(i, samp_af_list, master_frame_AF, master_frame_gen1, master_frame_gen2, pop_matrix, snpsamps_AF)
    master_frame_final
    
}
    # CHANGE WHERE YOU WANT TO WRITE FILE
    # TODO: generalize the name of file to include parameter variables
    #write.table(master_frame_final, 'mixtures/testfile2.txt', row.names = FALSE, quote = FALSE, sep = '\t')
