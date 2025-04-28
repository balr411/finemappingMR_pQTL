#This script extends the example in 00_02.PCSK9_PLA2G7_trans_analysis.R to use trans signals
#and cis signals for the PCSK9 protein

library("data.table")
library("stringr")
library("finemappingMR")
library("reticulate")
library("dplyr")
library("GenomicRanges")
library("Repitools")

use_condaenv("r-reticulate", required = TRUE) #Not sure if this will work great on the cluster?


df_res <- fread(cmd = "zcat /net/1000g/hmkang/etc/ukb/ppp/finemapping/sentinel/EURv2_merged_sentinel.cis_trans.wgs.tsv.gz",
                header = TRUE,
                data.table = FALSE)

#Filter to only PCSK9
df_res <- df_res[df_res$GENE == "PCSK9",]

#Remove HLA variants 
idx_hla <- which(df_res[[1]] == 6 & ((df_res$POS >= 28510120 & df_res$POS <= 33480577)))
df_res[idx_hla,] #There is one SNP at 6:31339002:G:A

if(length(idx_hla) > 0){
  df_res <- df_res[-idx_hla, ]
}

#Now remove variants with p-value < 5e-8
df_res <- df_res[df_res$LOG10P > -log10(5e-8),]

#Now create regions +- 500kb of these signals:
df_res$START <- pmax(df_res$POS - 500000, 1)
df_res$END <- df_res$POS + 500000

myranges <- GRanges(seqnames = df_res[[1]],
                    ranges = IRanges(start = df_res$START,
                                     end = df_res$END))

myranges_red <- reduce(myranges)
myranges_red_df <- annoGR2DF(myranges_red)

dim(myranges_red_df) #17 4 - Note that this is the same as the number of IVs (I think Hyun created them this way)

##################################################
#Now get a list of the ranges of the LD files that Hyun has at /net/1000g/hmkang/data/UKBB_LD/

#Note that some of the LD files end in the prefix ".npz2"
#According to the Alkes Price group (https://github.com/omerwe/polyfun/issues/17),
#These files contain regions where there was extremely long LD, and both SuSiE and 
#FINEMAP reported many false positives. Potentially it might be a good idea to delete
#them?

ld_files_gz <- system("ls /net/1000g/hmkang/data/UKBB_LD/", intern = TRUE)
ld_files_gz <- ld_files_gz[!(ld_files_gz %in% c("baselineLF_v2.2.UKB.polyfun.tar.gz", "baselineLF_v2.2.UKB.tar.gz", "readme_ld.txt"))]
ld_files_gz <- ld_files_gz[grepl("\\.gz$", ld_files_gz)]
ld_files_gz <- gsub(".gz", "", ld_files_gz)

tmp_vec <- unlist(strsplit(ld_files_gz, split = "_", fixed = TRUE))
chr_ld <- as.numeric(gsub("chr", "", tmp_vec[c(TRUE, FALSE, FALSE)]))
start_ld <- as.numeric(tmp_vec[c(FALSE, TRUE, FALSE)])
end_ld <- as.numeric(tmp_vec[c(FALSE, FALSE, TRUE)])

#Read in the summary statistics first:
source_python("../../Scripts/01.initial_analysis/00_01.python_scripts.py")

pcsk9_xty <- list()
pla2g7_xty <- list()
pcsk9_xtx <- list()
pla2g7_xtx <- list()

#Count the .npz2 files
npz2_indices <- c()

for(i in 1:nrow(myranges_red_df)){
  print(i)
  
  chr_curr <- as.integer(as.character(myranges_red_df$chr))[i] #For some reason this was a factor and was converting 15 to 9?? - should check that this doesn't happen in the other analysis
  reg_start <- as.integer(as.character(myranges_red_df$start))[i]
  reg_end <- as.integer(as.character(myranges_red_df$end))[i]
  
  df_full_sum_stat <- fread(cmd = str_glue("tabix -h /net/1000g/hmkang/etc/ukb/ppp/data/merged/EURv2.merged.{chr_curr}.tsv.gz {chr_curr}:{reg_start}-{reg_end}"),
                            header = TRUE,
                            data.table = FALSE)
  
  genes_full <- unlist(strsplit(colnames(df_full_sum_stat)[seq(from = 11, to = ncol(df_full_sum_stat), by = 6)], split = ":", fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)]
  idx_pcsk9 <- which(genes_full == "PCSK9") 
  idx_pla2g7 <- which(genes_full == "PLA2G7") 
  
  df_pcsk9 <- df_full_sum_stat[,c(1:5, (idx_pcsk9*6):(idx_pcsk9*6 + 5))]
  df_pla2g7 <- df_full_sum_stat[,c(1:5, (idx_pla2g7*6):(idx_pla2g7*6 + 5))]
  
  #Remove duplicated positions
  df_pcsk9 <- anti_join(df_pcsk9, df_pcsk9[duplicated(df_pcsk9[,2]),], by = "POS")
  df_pla2g7 <- anti_join(df_pla2g7, df_pla2g7[duplicated(df_pla2g7[,2]),], by = "POS")
  
  #Get rid of the rows with no beta:
  df_pcsk9 <- df_pcsk9[complete.cases(df_pcsk9),]
  df_pla2g7 <- df_pla2g7[complete.cases(df_pla2g7),]
  
  #Note for i = 1:
  #min(df_pcsk9[[6]]) 0.000733452 - so there are some pretty rare variants
  
  ############################################################################################
  #Remove variants with MAF < 0.01:
  maf_pcsk9 <- pmin(df_pcsk9[,6], 1 - df_pcsk9[,6])
  maf_pla2g7 <- pmin(df_pla2g7[,6], 1 - df_pla2g7[,6])
  
  df_pcsk9 <- df_pcsk9[maf_pcsk9 > 0.01,] #### Might not be correct !
  df_pla2g7 <- df_pla2g7[maf_pla2g7 > 0.01,] #### Might not be correct !
  ############################################################################################
  
  snps_common <- intersect(df_pcsk9$ID, df_pla2g7$ID)
  
  df_pcsk9 <- df_pcsk9[df_pcsk9$ID %in% snps_common,]
  df_pla2g7 <- df_pla2g7[df_pla2g7$ID %in% snps_common,]
  
  if(!isTRUE(all.equal(df_pcsk9$ID, df_pla2g7$ID))){
    print("Error!")
    break
  }
  
  #For i = 8, remove the one SNP that is lifted over weirdly (was i = 7 in trans analysis so should be i = 8 here)
  if(i == 8){
    df_pcsk9 <- df_pcsk9[-301,]
    df_pla2g7 <- df_pla2g7[-301,]
  }
  
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
  pcsk9_xty[[i]] <- sqrt(pcsk9_n - 1)*df_pcsk9$z_adj
  pla2g7_xty[[i]] <- sqrt(pla2g7_n - 1)*df_pla2g7$z_adj
  
  #Now try to match the LD
  pos_b37 <- as.numeric(unlist(strsplit(df_pcsk9$ID, split = ":", fixed = TRUE))[c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)])
  lower_b37 <- min(pos_b37)
  upper_b37 <- max(pos_b37)
  
  idx_ld <- which(chr_ld == chr_curr & lower_b37 >= start_ld & upper_b37 <= end_ld)
  
  if(length(idx_ld) == 0){
    print(str_glue("LD is mismatched for i = {i}"))
  }
  
  #Now read in the LD
  #Note that there could be multiple LD files containing the region, so just choose 
  #the first one
  reg_start_b37 <- start_ld[idx_ld[1]]
  reg_end_b37 <- end_ld[idx_ld[1]]
  
  ld_pos <- fread(cmd = str_glue("zcat /net/1000g/hmkang/data/UKBB_LD/chr{chr_curr}_{reg_start_b37}_{reg_end_b37}.gz"),
                  header = TRUE,
                  data.table = FALSE)
  
  #Some positions may not appear in the LD matrix, so check that first:
  pcsk9_id_rd <- gsub(":imp:v1", "", df_pcsk9$ID, fixed = TRUE)
  ld_pos_id <- paste(ld_pos$chromosome, ld_pos$position, ld_pos$allele1, ld_pos$allele2, sep = ":")
  
  common_snps_ld <- intersect(pcsk9_id_rd, ld_pos_id)
  
  #Now get indices:
  idx_to_output_ld <- which(ld_pos_id %in% pcsk9_id_rd) ## Duplicate positions (multi-allelic variants) will mess this up - go back and delete them (not really sure why this is true?)
  
  if(!isTRUE(all.equal(ld_pos_id[idx_to_output_ld], pcsk9_id_rd))){
    print("LD variants don't match PCSK9 variants!")
    idx_common_sumstat <- which(pcsk9_id_rd %in% common_snps_ld)
    pcsk9_xty[[i]] <- pcsk9_xty[[i]][idx_common_sumstat]
    pla2g7_xty[[i]] <- pla2g7_xty[[i]][idx_common_sumstat]
    idx_to_output_ld <- which(ld_pos_id %in% common_snps_ld)
    
    if(!isTRUE(all.equal(ld_pos_id[idx_to_output_ld], pcsk9_id_rd[idx_common_sumstat]))){
      print("LD variants still don't match PCSK9 variants!")
      break
    }
  }
  
  
  #Python indices start at 0 so subtract one from each index: 
  idx_to_output_ld <- idx_to_output_ld - 1
  ld_name_curr <- str_glue("/net/1000g/hmkang/data/UKBB_LD/chr{chr_curr}_{reg_start_b37}_{reg_end_b37}.npz")
  
  #Check to see if the .npz or .npz2 files exist
  if(!file.exists(ld_name_curr)){
    ld_name_curr <- str_glue("/net/1000g/hmkang/data/UKBB_LD/chr{chr_curr}_{reg_start_b37}_{reg_end_b37}.npz2")
    npz2_indices <- c(npz2_indices, i)
  }
  
  
  ld_mat_red <- sparse_ld_subset(ld_name_curr, idx_to_output_ld)
  
  
  pcsk9_xtx[[i]] <- (pcsk9_n - 1)*ld_mat_red
  pla2g7_xtx[[i]] <- (pla2g7_n - 1)*ld_mat_red
  
  
}

npz2_indices

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




