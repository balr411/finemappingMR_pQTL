#This script will just try to read in and perform some analyses with Hyun's files
#Can first start by trying to read in the sentinel variants and getting a region 
#ready to run 

library("data.table")
library("stringr")
#library("RcppCNPy")
library("finemappingMR")
library("reticulate")
library("dplyr")

use_condaenv("r-reticulate", required = TRUE) #Not sure if this will work great on the cluster?


df_res <- fread("/net/1000g/hmkang/etc/ukb/ppp/finemapping/sentinel/EURv2_merged_sentinel.cis_trans.wgs.tsv.gz",
                header = TRUE,
                data.table = FALSE)

#Start with PCSK9 
#The cis signal is 1:55505647:G:T:imp:v1 in build 37
#In build 38, it is 1:55039974:G:T

#Note that Hyun then uses the /net/1000g/hmkang/etc/ukb/ppp/data/merged/EURv2.merged.1.tsv.gz
#file in order to run ACAT to check if there are significant trans signals
#for the variant 1:55505647:G:T:imp:v1

#The script /net/1000g/hmkang/etc/ukb/ppp/scripts/test_trans_inflation.py gathers 
#the p-value for the variant 1:55039974:G:T (or just simply 1:55039974) across 
#all of the other proteins

#Read in the results for this particular variant:

df_pcsk9_cis_var <- fread(cmd = "tabix -h /net/1000g/hmkang/etc/ukb/ppp/data/merged/EURv2.merged.1.tsv.gz 1:55039974-55039974",
                          header = TRUE,
                          data.table = FALSE)

num_genes <- (ncol(df_pcsk9_cis_var) - 5)/6

num_genes #2940

#See which variant has the highest signal and use this one first
genes <- unlist(strsplit(colnames(df_pcsk9_cis_var)[seq(from = 11, to = ncol(df_pcsk9_cis_var), by = 6)], split = ":", fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)]
genes_logp <- as.numeric(unname(df_pcsk9_cis_var[1, seq(from = 11, to = ncol(df_pcsk9_cis_var), by = 6)]))

genes_logp <- genes_logp[genes != "PCSK9"]
genes <- genes[genes != "PCSK9"]

which.max(genes_logp) #2075
genes[2075] #PLA2G7
genes_logp[2075] #15.2843

#Note that PLA2G7 binds to LDL so this makes sense
#PCSK9 has a positive effect on LDL, and the variant 1:55039974:G:T has a negative
#beta for both PCSK9 and PLA2G7, which makes sense.

#Start by running finemappingMR on the region +- 500kb of 1:55039974:G:T on PCSK9 -> PLA2G7
chr_curr <- 1
reg_start <- 55039974 - 500000
reg_end <- 55039974 + 500000

#Now read in the summary statistics
df_full_sum_stat <- fread(cmd = str_glue("tabix -h /net/1000g/hmkang/etc/ukb/ppp/data/merged/EURv2.merged.{chr_curr}.tsv.gz {chr_curr}:{reg_start}-{reg_end}"),
                          header = TRUE,
                          data.table = FALSE)

dim(df_full_sum_stat)
#6680 17645

## I think this is actually not the correct file. This one probably only has variants 
#that were significant for at least one protein?

#Double check by reading in the files for PCSK9 and PLA2G7 separately:
df_pcsk9 <- fread(cmd = str_glue("tabix -h /net/1000g/hmkang/data/UKB/PPP/sumstats/reformatted/EURv2/EURv2_PCSK9_Q8NBP7_OID20235_v1_Cardiometabolic.gz {chr_curr}:{reg_start}-{reg_end}"),
                  header = TRUE,
                  data.table = FALSE)

dim(df_pcsk9)
#6665   14

df_pla2g7 <- fread(cmd = str_glue("tabix -h /net/1000g/hmkang/data/UKB/PPP/sumstats/reformatted/EURv2/EURv2_PLA2G7_Q13093_OID21105_v1_Neurology.gz {chr_curr}:{reg_start}-{reg_end}"),
                   header = TRUE,
                   data.table = FALSE)

dim(df_pla2g7)
#6658   14

length(intersect(df_pla2g7$ID, df_full_sum_stat$ID)) #6658
length(intersect(df_pcsk9$ID, df_full_sum_stat$ID)) #6665

#Hence it does indeed contain all of the variants. Proceed to getting the summary stats:
genes_full <- unlist(strsplit(colnames(df_full_sum_stat)[seq(from = 11, to = ncol(df_full_sum_stat), by = 6)], split = ":", fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)]
which(genes_full == "PCSK9") #2004
which(genes_full == "PLA2G7") #2076

df_pcsk9 <- df_full_sum_stat[,c(1:5, (2004*6):(2004*6 + 5))]
df_pla2g7 <- df_full_sum_stat[,c(1:5, (2076*6):(2076*6 + 5))]

#Remove duplicated positions
df_pcsk9 <- anti_join(df_pcsk9, df_pcsk9[duplicated(df_pcsk9[,2]),], by = "POS")
df_pla2g7 <- anti_join(df_pla2g7, df_pla2g7[duplicated(df_pla2g7[,2]),], by = "POS")


#Get rid of the rows with no beta:
df_pcsk9 <- df_pcsk9[complete.cases(df_pcsk9),]
df_pla2g7 <- df_pla2g7[complete.cases(df_pla2g7),]

#Remove variants with MAF < 0.01:
df_pcsk9 <- df_pcsk9[df_pcsk9[,6] > 0.01,]
df_pla2g7 <- df_pla2g7[df_pla2g7[,6] > 0.01,]

snps_common <- intersect(df_pcsk9$ID, df_pla2g7$ID)

df_pcsk9 <- df_pcsk9[df_pcsk9$ID %in% snps_common,]
df_pla2g7 <- df_pla2g7[df_pla2g7$ID %in% snps_common,]

nrow(df_pcsk9) #4113 (was 4124 before removing multi-allelics)
nrow(df_pla2g7) #4113

all.equal(df_pcsk9$ID, df_pla2g7$ID) #TRUE

#Now adjust the z-scores as in the SuSiE code (i.e. PVE-adjusted Z-scores):
pcsk9_n <- as.numeric(df_pcsk9[1,7])
pcsk9_adj <- (pcsk9_n - 1)/(df_pcsk9[,10]^2 + pcsk9_n - 2)
z_pcsk9   <- sqrt(pcsk9_adj) * df_pcsk9[,10]
df_pcsk9$z_adj <- z_pcsk9

pla2g7_n <- as.numeric(df_pla2g7[1,7])
pla2g7_adj <- (pla2g7_n - 1)/(df_pla2g7[,10]^2 + pla2g7_n - 2)
z_pla2g7   <- sqrt(pla2g7_adj) * df_pla2g7[,10]
df_pla2g7$z_adj <- z_pla2g7

#Results given on standardized scale if var(y) unknown(?) - again taken from SuSiE code
pcsk9_xty <- list()
pcsk9_xty[[1]] <- sqrt(pcsk9_n - 1)*df_pcsk9$z_adj

pla2g7_xty <- list()
pla2g7_xty[[1]] <- sqrt(pla2g7_n - 1)*df_pla2g7$z_adj

##################################################################################
#Now read in the LD
#First convert the ranges to build 37
pos_b37 <- as.numeric(unlist(strsplit(df_pcsk9$ID, split = ":", fixed = TRUE))[c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)])
lower_b37 <- min(pos_b37)
upper_b37 <- max(pos_b37)

#Region is 55005704:56005610
#Hence should be able to use the file /net/1000g/hmkang/data/UKBB_LD/chr1_55000001_58000001.npz

#To read in the LD, I am first going to read in the file /net/1000g/hmkang/data/UKBB_LD/chr1_55000001_58000001.gz,
#and match the indices that I'll need LD for to the indices in that file. Then
#I will output these indices, and use them to write out a new (much smaller) .npz
#file which can be converted to dense format in R

ld_pos <- fread("/net/1000g/hmkang/data/UKBB_LD/chr1_55000001_58000001.gz",
                header = TRUE,
                data.table = FALSE)

#Some positions may not appear in the LD matrix, so check that first:
pcsk9_id_rd <- gsub(":imp:v1", "", df_pcsk9$ID, fixed = TRUE)
ld_pos_id <- paste(ld_pos$chromosome, ld_pos$position, ld_pos$allele1, ld_pos$allele2, sep = ":")
sum(ld_pos_id %in% pcsk9_id_rd) #4113 hence they are all present

#Now get indices:
idx_to_output <- which(ld_pos_id %in% pcsk9_id_rd) ## Duplicate positions (multi-allelic variants) will mess this up - go back and delete them
all.equal(ld_pos_id[idx_to_output], pcsk9_id_rd) #TRUE

#np <- import("numpy") -> This causes issues when calling sparse_ld_subset(), because the ndarray is not automatically converted to a matrix
#In the future if I want to use something like this, I need to add convert = FALSE to the import() command
#ld_mat <- np$load("/net/1000g/hmkang/data/UKBB_LD/chr1_55000001_58000001.npz")
source_python("../../Scripts/01.initial_analysis/00_01.python_scripts.py")

#Python indices start at 0 so subtract one from each index: 
idx_to_output <- idx_to_output - 1
ld_mat_red <- sparse_ld_subset("/net/1000g/hmkang/data/UKBB_LD/chr1_55000001_58000001.npz", idx_to_output)

## This works

pcsk9_xtx <- list()
pcsk9_xtx[[1]] <- (pcsk9_n - 1)*ld_mat_red

pla2g7_xtx <- list()
pla2g7_xtx[[1]] <- (pla2g7_n - 1)*ld_mat_red

pcsk9_var <- 1
pla2g7_var <- 1

#Now call finemappingMR 
res_our_method <- run_freq_method_ss(pcsk9_xtx, pcsk9_xty, pcsk9_n - 1,
                                     pla2g7_xtx, pla2g7_xty, pla2g7_n - 1,
                                     pcsk9_n, pla2g7_n,
                                     L_x = 10, L_y = 10,
                                     scaled_prior_variance_x = 0.2,
                                     scaled_prior_variance_y = 0.2,
                                     estimate_prior_variance_x = TRUE,
                                     estimate_prior_variance_y = TRUE,
                                     residual_variance_x = NULL,
                                     residual_variance_y = NULL,
                                     estimate_residual_variance_x = FALSE,
                                     estimate_residual_variance_y = FALSE,
                                     tol = 1e-4,
                                     max_iter = 1000,
                                     calc_cs_x = TRUE,
                                     calc_cs_y = TRUE,
                                     verbose = TRUE)



res_our_method$res
#gamma    gamma_var iter
#0.1991038 0.0004407491   14

2*pnorm(-0.1991038/sqrt(0.0004407491)) #2.451137e-21

#Highly significant and in the correct direction


