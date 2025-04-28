#This script looks at the PCSK9-PLA2G7 example using only trans signals, but runs
#finemappingMR one at a time

library("data.table")
library("stringr")
library("finemappingMR")
library("reticulate")
library("dplyr")
library("GenomicRanges")
library("Repitools")
library("ggplot2")
library("TwoSampleMR")

#Modified RAPS function 
mr_modified <- function (dat, 
                         parameters = default_parameters(), 
                         method_list = subset(mr_method_list(), use_by_default)$obj){
  
  mr_raps_modified <- function (b_exp, b_out, se_exp, se_out,parameters) 
  {
    out <- try(suppressMessages(mr.raps::mr.raps(b_exp, b_out, se_exp, se_out,
                                                 over.dispersion = parameters$over.dispersion, 
                                                 loss.function = parameters$loss.function,
                                                 diagnosis = FALSE)),
               silent = T)
    
    # The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion
    # When encountering such warning, change the over.dispersion as 'FASLE'
    
    if ('try-error' %in% class(out))
    {
      output = list(b = NA, se = NA, pval = NA, nsnp = NA)
    }
    else
    {
      output = list(b = out$beta.hat, se = out$beta.se, 
                    pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 2, nsnp = length(b_exp))
    }
    return(output)
  }
  
  method_list_modified <- stringr::str_replace_all(method_list, "mr_raps","mr_raps_modified")
  
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"),function(x1)
  {
    x <- subset(x1, mr_keep)
    
    if (nrow(x) == 0) {
      message("No SNPs available for MR analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    }
    else {
      message("Analysing '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
    }
    res <- lapply(method_list_modified, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    }
    )
    
    methl <- mr_method_list()
    mr_tab <- data.frame(outcome = x$outcome[1], exposure = x$exposure[1], 
                         method = methl$name[match(method_list, methl$obj)], 
                         nsnp = sapply(res, function(x) x$nsnp), 
                         b = sapply(res, function(x) x$b), 
                         se = sapply(res, function(x) x$se), 
                         pval = sapply(res, function(x) x$pval))
    
    mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
    
    return(mr_tab)
  }
  )
  return(mr_tab)
}


use_condaenv("r-reticulate", required = TRUE) #Not sure if this will work great on the cluster?


df_res <- fread(cmd = "zcat /net/1000g/hmkang/etc/ukb/ppp/finemapping/sentinel/EURv2_merged_sentinel.cis_trans.wgs.tsv.gz",
                header = TRUE,
                data.table = FALSE)

#Filter to only PCSK9
df_res <- df_res[df_res$GENE == "PCSK9",]

#Filter to only include trans signals: 
df_res_trans <- df_res[df_res$TYPE == "trans",]

#Remove HLA variants 
idx_hla <- which(df_res_trans[[1]] == 6 & ((df_res_trans$POS >= 28510120 & df_res_trans$POS <= 33480577)))
df_res_trans[idx_hla,] #There is one SNP at 6:31339002:G:A

if(length(idx_hla) > 0){
  df_res_trans <- df_res_trans[-idx_hla, ]
}

#Now remove variants with p-value < 5e-8
df_res_trans <- df_res_trans[df_res_trans$LOG10P > -log10(5e-8),]

#Now create regions +- 500kb of these signals:
df_res_trans$START <- pmax(df_res_trans$POS - 500000, 1)
df_res_trans$END <- df_res_trans$POS + 500000

myranges <- GRanges(seqnames = df_res_trans[[1]],
                    ranges = IRanges(start = df_res_trans$START,
                                     end = df_res_trans$END))

myranges_red <- reduce(myranges)
myranges_red_df <- annoGR2DF(myranges_red)

dim(myranges_red_df) #16 4 - Note that this is the same as the number of IVs (I think Hyun created them this way)

##################################################
#Now get a list of the ranges of the LD files that Hyun has at /net/1000g/hmkang/data/UKBB_LD/
#ld_files <- system("ls /net/1000g/hmkang/data/UKBB_LD/", intern = TRUE)
#ld_files <- ld_files[!(ld_files %in% c("baselineLF_v2.2.UKB.polyfun.tar.gz", "baselineLF_v2.2.UKB.tar.gz", "readme_ld.txt"))]
#ld_files <- ld_files[!grepl("\\.gz$", ld_files)]
#ld_files <- gsub(".npz", "", ld_files)

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


#Note that this is quite a lot of regions so I am not sure how much memory this will require
#Read in the summary statistics first:
source_python("../../Scripts/01.initial_analysis/00_01.python_scripts.py")

pcsk9_xty <- list()
pla2g7_xty <- list()
pcsk9_xtx <- list()
pla2g7_xtx <- list()

res_our_method_list <- list()

#df_pla2g7_sent <- data.frame()

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
  
  #For i = 7, remove the one SNP that is lifted over weirdly
  if(i == 7){
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
  
  #Keep track of the sentinel variants of pcsk9 in pla2g7 as well
  #df_pla2g7_sent <- rbind(df_pla2g7_sent, df_pla2g7[df_pla2g7$ID == df_res_trans$SNPID[i],])
  
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
  
  
  pcsk9_xtx[[i]] <- (pcsk9_n - 1)*ld_mat_red #could n differ across regions?
  pla2g7_xtx[[i]] <- (pla2g7_n - 1)*ld_mat_red
  
  pcsk9_var <- 1
  pla2g7_var <- 1
  
  #Now call finemappingMR 
  pcsk9_xtx_list_curr <- list()
  pcsk9_xtx_list_curr[[1]] <- pcsk9_xtx[[i]]
  
  pcsk9_xty_list_curr <- list()
  pcsk9_xty_list_curr[[1]] <- pcsk9_xty[[i]]
  
  pla2g7_xtx_list_curr <- list()
  pla2g7_xtx_list_curr[[1]] <- pla2g7_xtx[[i]]
  
  pla2g7_xty_list_curr <- list()
  pla2g7_xty_list_curr[[1]] <- pla2g7_xty[[i]]
  
  res_our_method <- run_freq_method_ss(pcsk9_xtx_list_curr, pcsk9_xty_list_curr, pcsk9_n - 1,
                                       pla2g7_xtx_list_curr, pla2g7_xty_list_curr, pla2g7_n - 1,
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
  
  res_our_method_list[[i]] <- res_our_method
  
}

npz2_indices #NULL

df_pla2g7_sent <- data.frame()
for(i in 1:nrow(df_res_trans)){
  
  print(i)
  
  chr_curr <- as.integer(as.character(myranges_red_df$chr))[i] #For some reason this was a factor and was converting 15 to 9?? - should check that this doesn't happen in the other analysis
  reg_start <- as.integer(as.character(myranges_red_df$start))[i]
  reg_end <- as.integer(as.character(myranges_red_df$end))[i]
  
  df_full_sum_stat <- fread(cmd = str_glue("tabix -h /net/1000g/hmkang/etc/ukb/ppp/data/merged/EURv2.merged.{chr_curr}.tsv.gz {chr_curr}:{reg_start}-{reg_end}"),
                            header = TRUE,
                            data.table = FALSE)
  
  genes_full <- unlist(strsplit(colnames(df_full_sum_stat)[seq(from = 11, to = ncol(df_full_sum_stat), by = 6)], split = ":", fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)]
  idx_pla2g7 <- which(genes_full == "PLA2G7") 
  df_pla2g7 <- df_full_sum_stat[,c(1:5, (idx_pla2g7*6):(idx_pla2g7*6 + 5))]
  
  df_pla2g7_sent <- rbind(df_pla2g7_sent, df_pla2g7[df_pla2g7$ID == df_res_trans$SNPID[i],])
  
  if(!any(df_pla2g7$ID == df_res_trans$SNPID[i])){
    print(i)
  }
  
}

#Get our gamma estimates and make a data frame
df_comb <- data.frame(pcsk9_beta = df_res_trans$BETA,
                      pcsk9_se = df_res_trans$SE,
                      pcsk9_Z = df_res_trans$Z,
                      pla2g7_beta = df_pla2g7_sent[[8]],
                      pla2g7_se = df_pla2g7_sent[[9]],
                      pla2g7_Z = df_pla2g7_sent[[10]],
                      fmr_gamma = NA,
                      fmr_se = NA,
                      fmr_Z = NA)

for(i in 1:length(res_our_method_list)){
  res_curr <- res_our_method_list[[i]]$res
  df_comb$fmr_gamma[i] <- res_curr$gamma[1]
  df_comb$fmr_se[i] <- sqrt(res_curr$gamma_var[1])
  df_comb$fmr_Z[i] <- df_comb$fmr_gamma[i]/df_comb$fmr_se[i]
}

#Add the other information needed for TwoSampleMR (go back later and double check that the alleles aren't swapped? - they appear not to be)
df_comb$effect_allele <- df_res_trans$ALT
df_comb$SNP <- df_res_trans$SNPID

#Put in the TwoSampleMR format
df_exposure <- data.frame(SNP = df_comb$SNP,
                          beta = df_comb$pcsk9_beta,
                          se = df_comb$pcsk9_se,
                          effect_allele = df_comb$effect_allele)

df_outcome <- data.frame(SNP = df_comb$SNP,
                         beta = df_comb$pla2g7_beta,
                         se = df_comb$pla2g7_se,
                         effect_allele = df_comb$effect_allele)

df_exposure <- format_data(df_exposure, type = "exposure")
df_outcome <- format_data(df_outcome, type = "outcome")

#Put action = 1 to make it ignore palindromic SNPs
dat_mr <- harmonise_data(
  exposure_dat = df_exposure, 
  outcome_dat = df_outcome,
  action = 1
)

#Note that orginally it tried to remove everything with the message:
#Removing the following SNPs for being palindromic with intermediate allele frequencies (listed all SNPs)

#Run traditional MR
df_raps_curr <- mr_modified(dat_mr, method_list = c("mr_raps"), 
                            parameters = list(over.dispersion = TRUE, loss.function = "tukey"))

#mr.raps::mr.raps.all(dat_mr$beta.exposure, dat_mr$beta.outcome, dat_mr$se.exposure, dat_mr$se.outcome) - gives different answers for huber compared to tukey

df_mr_curr <- mr(dat_mr, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

#Run SMR
pval_exp <- dat_mr$pval.exposure
max_snp <- which.min(pval_exp)
Z_x <- dat_mr$beta.exposure[max_snp]/dat_mr$se.exposure[max_snp]
Z_y <- dat_mr$beta.outcome[max_snp]/dat_mr$se.outcome[max_snp]

test_stat <- (Z_x^2 * Z_y^2)/(Z_x^2 + Z_y^2)
p_smr <- pchisq(test_stat, 1, lower.tail = FALSE)

smr_est <- dat_mr$beta.outcome[max_snp]/dat_mr$beta.exposure[max_snp]
smr_est_var <- ((dat_mr$beta.outcome[max_snp]^2)/(dat_mr$beta.exposure[max_snp]^2)) * (dat_mr$se.outcome[max_snp]^2/dat_mr$beta.outcome[max_snp]^2 + dat_mr$se.exposure[max_snp]^2/dat_mr$beta.exposure[max_snp]^2)



## Check inverse-variance weighted estimate of the regions
gamma_ests <- df_comb$fmr_gamma
gamma_vars <- df_comb$fmr_se^2

wi <- 1/gamma_vars
se_curr <- sqrt(1/sum(wi))
gamma_comb <- sum(gamma_ests*wi)/sum(wi)

Z_curr <- gamma_comb/se_curr

## Now make a MR plot of the instruments
#Now add the error bars and make the plot
df_comb$Y_err <- df_comb$pla2g7_se * 1.96
df_comb$X_err <- df_comb$pcsk9_se * 1.96

#Find limits to use in the plot
max(df_comb$pcsk9_beta + df_comb$X_err) #0.3611869
min(df_comb$pcsk9_beta - df_comb$X_err) #-0.1522346

max(df_comb$pla2g7_beta + df_comb$Y_err) #0.06455847
min(df_comb$pla2g7_beta - df_comb$Y_err) #-0.5061643

cols <- c("RAPS" = "#1B9E77", "IVW"= "#D95F02", "Median"= "#7570B3", "SMR" = "#E7298A", "finemappingMR" = "#66A61E")
maf_colors <- c("<0.05" = "blue", ">0.05" = "red")

p1 <- ggplot(df_comb, aes(x = pcsk9_beta, y = pla2g7_beta)) + geom_point() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.position = "right",
        axis.line = element_line(color='black'),
        plot.margin = margin(0,0,0,0, 'cm')) +  
  xlim(-0.51, 0.51) + ylim(-0.51, 0.51) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_abline(aes(intercept = 0, slope = -0.037, color = "finemappingMR"))  + 
  geom_abline(aes(intercept = 0, slope = -0.1951792, color = "SMR")) +
  geom_abline(aes(intercept = 0, slope = -0.06344671, color = "RAPS")) +
  geom_abline(aes(intercept = 0, slope = 0.3779086, color = "IVW")) +
  geom_abline(aes(intercept = 0, slope = -0.1075358, color = "Median")) +
  ggtitle("Scatter Plot of PLA2G7 Betas vs. PCSK9 Betas Using Trans Signals Only") + 
  ylab("PLA2G7 Betas") + xlab("PCSK9 Betas") + 
  geom_errorbar(aes(ymin = pla2g7_beta - Y_err, ymax = pla2g7_beta + Y_err), width = 0) + 
  geom_errorbarh(aes(xmin = pcsk9_beta - X_err, xmax = pcsk9_beta + X_err), height = 0) +
  scale_color_manual(name = "Method", values = cols, guide = "legend") 

png("MR_plot_PLA2G7_PCSK9_trans_only.png", width = 600, height = 600)
p1
dev.off()

###################################################################################
#Now do the cis + trans analysis for the traditional methods
chr_curr <- 1
reg_start <- 55039970
reg_end <- 55039980

df_full_sum_stat <- fread(cmd = str_glue("tabix -h /net/1000g/hmkang/etc/ukb/ppp/data/merged/EURv2.merged.{chr_curr}.tsv.gz {chr_curr}:{reg_start}-{reg_end}"),
                          header = TRUE,
                          data.table = FALSE)

genes_full <- unlist(strsplit(colnames(df_full_sum_stat)[seq(from = 11, to = ncol(df_full_sum_stat), by = 6)], split = ":", fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)]
idx_pla2g7 <- which(genes_full == "PLA2G7") 
df_pla2g7 <- df_full_sum_stat[,c(1:5, (idx_pla2g7*6):(idx_pla2g7*6 + 5))]

df_pla2g7_sent_cis_trans <- rbind(df_pla2g7_sent, df_pla2g7[df_pla2g7$ID == "1:55505647:G:T:imp:v1",])

#Change the MR data frames
#Put in the TwoSampleMR format
df_exposure <- data.frame(SNP = df_comb$SNP,
                          beta = df_comb$pcsk9_beta,
                          se = df_comb$pcsk9_se,
                          effect_allele = df_comb$effect_allele)

df_outcome <- data.frame(SNP = df_comb$SNP,
                         beta = df_comb$pla2g7_beta,
                         se = df_comb$pla2g7_se,
                         effect_allele = df_comb$effect_allele)

#Add the cis-region 
df_exposure <- rbind(df_exposure, data.frame(SNP = df_res$SNPID[2],
                                             beta = df_res$BETA[2],
                                             se = df_res$SE[2],
                                             effect_allele = df_res$ALT[2]))

df_outcome <- rbind(df_outcome, data.frame(SNP = df_res$SNPID[2],
                                           beta = df_pla2g7_sent_cis_trans[17,8],
                                           se = df_pla2g7_sent_cis_trans[17,9],
                                           effect_allele = df_res$ALT[2]))

df_exposure <- format_data(df_exposure, type = "exposure")
df_outcome <- format_data(df_outcome, type = "outcome")

#Put action = 1 to make it ignore palindromic SNPs
dat_mr <- harmonise_data(
  exposure_dat = df_exposure, 
  outcome_dat = df_outcome,
  action = 1
)

#Now get the SMR estimate for just the cis-region (note this is the strongest signal anyways so it is also the SMR estimate for cis + trans)
max_snp <- 1
Z_x <- dat_mr$beta.exposure[max_snp]/dat_mr$se.exposure[max_snp]
Z_y <- dat_mr$beta.outcome[max_snp]/dat_mr$se.outcome[max_snp]

test_stat <- (Z_x^2 * Z_y^2)/(Z_x^2 + Z_y^2)
p_smr <- pchisq(test_stat, 1, lower.tail = FALSE)

smr_est <- dat_mr$beta.outcome[max_snp]/dat_mr$beta.exposure[max_snp]
smr_est_var <- ((dat_mr$beta.outcome[max_snp]^2)/(dat_mr$beta.exposure[max_snp]^2)) * (dat_mr$se.outcome[max_snp]^2/dat_mr$beta.outcome[max_snp]^2 + dat_mr$se.exposure[max_snp]^2/dat_mr$beta.exposure[max_snp]^2)

#Run traditional MR
df_raps_curr <- mr_modified(dat_mr, method_list = c("mr_raps"), 
                            parameters = list(over.dispersion = TRUE, loss.function = "tukey"))

#mr.raps::mr.raps.all(dat_mr$beta.exposure, dat_mr$beta.outcome, dat_mr$se.exposure, dat_mr$se.outcome) - gives different answers for huber compared to tukey

df_mr_curr <- mr(dat_mr, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

## Now make a MR plot of the instruments
#Add the cis- snp (ignore all of the other parts and just change beta, se)
df_comb <- rbind(df_comb, df_comb[16,])

df_comb$pcsk9_beta[17] <- dat_mr$beta.exposure[1]
df_comb$pcsk9_se[17] <- dat_mr$se.exposure[1]
df_comb$pla2g7_beta[17] <- dat_mr$beta.outcome[1]
df_comb$pla2g7_se[17] <- dat_mr$se.outcome[1]

#Now add the error bars and make the plot
df_comb$Y_err <- df_comb$pla2g7_se * 1.96
df_comb$X_err <- df_comb$pcsk9_se * 1.96

#Find limits to use in the plot
max(df_comb$pcsk9_beta + df_comb$X_err) #0.3611869
min(df_comb$pcsk9_beta - df_comb$X_err) #-1.20774

max(df_comb$pla2g7_beta + df_comb$Y_err) #0.06455847
min(df_comb$pla2g7_beta - df_comb$Y_err) #-0.5061643

cols <- c("RAPS" = "#1B9E77", "IVW"= "#D95F02", "Median"= "#7570B3", "SMR" = "#E7298A", "finemappingMR" = "#66A61E")
maf_colors <- c("<0.05" = "blue", ">0.05" = "red")

p1 <- ggplot(df_comb, aes(x = pcsk9_beta, y = pla2g7_beta)) + geom_point() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        plot.title = element_text(size = 16, hjust = 0.5), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.position = "right",
        axis.line = element_line(color='black'),
        plot.margin = margin(0,0,0,0, 'cm')) +  
  xlim(-1.21, 1.21) + ylim(-1.21, 1.21) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_abline(aes(intercept = 0, slope = 0.159075, color = "finemappingMR"))  + 
  geom_abline(aes(intercept = 0, slope = 0.1901732, color = "SMR")) +
  geom_abline(aes(intercept = 0, slope = 0.1754628, color = "RAPS")) +
  geom_abline(aes(intercept = 0, slope = 0.2522063, color = "IVW")) +
  geom_abline(aes(intercept = 0, slope = 0.1778820, color = "Median")) +
  ggtitle("Scatter Plot of PLA2G7 Betas vs. PCSK9 Betas Using Cis and Trans") + 
  ylab("PLA2G7 Betas") + xlab("PCSK9 Betas") + 
  geom_errorbar(aes(ymin = pla2g7_beta - Y_err, ymax = pla2g7_beta + Y_err), width = 0) + 
  geom_errorbarh(aes(xmin = pcsk9_beta - X_err, xmax = pcsk9_beta + X_err), height = 0) +
  scale_color_manual(name = "Method", values = cols, guide = "legend") 

png("MR_plot_PLA2G7_PCSK9_cis_and_trans.png", width = 600, height = 600)
p1
dev.off()


