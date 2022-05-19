library(dplyr)
library(tidyr)
library(stringr)


#read table
p53.gnomAD <- read.csv("/Users/elpap/Dropbox (ELELAB)/ELELAB Team Folder/2021.WORK/PAPER_IN_PREPARATION/p53_DBD_CANCER_MUTATIONS/data/gnomAD/p53_gnomad.csv", header = TRUE,sep =",")
#332 
p53.gnomAD.missense <- subset(p53.gnomAD, major_consequence=="missense_variant")
#217
p53.gnomAD.missense$protein_position <- as.numeric(p53.gnomAD.missense$protein_position)
#select only DBD
p53.gnomAD.missense.DBD <- subset(p53.gnomAD.missense, ((protein_position >=98) & (protein_position <=289)))
#116
write.table(p53.gnomAD.missense.DBD, file="/Users/elpap/Dropbox (ELELAB)/ELELAB Team Folder/2021.WORK/PAPER_IN_PREPARATION/p53_DBD_CANCER_MUTATIONS/data/gnomAD/p53_missense_DBD_gnomad.txt", quote=FALSE)
mut_list_p53_gnomAD_DBD <-select(p53.gnomAD.missense.DBD,wt_residue,protein_position,mutated_residue)
write.table(mut_list_p53_gnomAD_DBD, file= "/Users/elpap/Dropbox (ELELAB)/ELELAB Team Folder/2021.WORK/PAPER_IN_PREPARATION/p53_DBD_CANCER_MUTATIONS/data/gnomAD/(mut_list_p53_gnomAD_DBD.txt", quote=FALSE)