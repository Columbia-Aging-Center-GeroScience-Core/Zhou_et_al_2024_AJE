# This document contains PGI computation code for Zhou et al. 2024 as of 05/14/2024
# Calculations followed the tutorial by Privé et al: https://privefl.github.io/bigsnpr/articles/LDpred2.html
#*Data sources*
# 1. DHWFS genomic data (target file)
# 2. Summary statistics for phenotype of interest (base file)    
# 3. HapMap3+ variants with independent LD blocks (https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061)

rm(list = ls())

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(tidyverse)
library(dplyr)
library(gee)
library(geepack)
library(data.table)
library(Matrix) # for sparse matrix

#### 1. Import Base file ------------------------------------------------------------------------------------------------------------------------------------------------
GWAS_sumstats = bigreadr::fread2("sumstats.txt.gz")
dim(GWAS_sumstats) 
GWAS_sumstats[1:5, ] 
sum(is.na(GWAS_sumstats))


### Load variants provided in the LD reference
# LD reference (correlations between pairs of genetic variants) for 1,054,330 HapMap3 variants based on 362,320 European individuals of the UK biobank.
# Pairs of variants further than 3 cM apart are assumed to have 0 correlation.
# REFERENCE: https://figshare.com/articles/dataset/European_LD_reference/13034123
map_ldref = readRDS("map.rds") # LD reference provided
dim(map_ldref) # 1,054,330 snps x 11 columns
map_ldref[1:5, ]
sum(is.na(map_ldref)) 

#### 2. Import Target file ------------------------------------------------------------------------------------------------------------------------------------------------
### Read from bed/bim/fam, it generates .bk and .rds files.
snp_readBed("HWFS_merged_tw.bed")
### Attach the "bigSNP" object in R session
obj_bigSNP = snp_attach("HWFS_merged_tw.rds")
map = setNames(obj_bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
G = obj_bigSNP$genotypes
NCORES = nb_cores()
obj_bigSNP$fam
str(obj_bigSNP, max.level = 2, strict.width = "cut")
dim(G) # 950 x 47,119,148
map[1:5, ]

### Match GWAS_sumstats and map_ldref
colnames(GWAS_sumstats) <- c('chr','rsid','pos', 'snp','a0', 'a1', 'beta','beta_se','n_eff', 'p', 'eaf')
GWAS_sumstats[1:5,]

df_beta1 = snp_match(GWAS_sumstats, map_ldref)

df_beta1[1:10, ]
df_beta1 = df_beta1[, c('chr', 'pos', 'a0', 'a1', 'beta', 'beta_se', 'n_eff', 'p', # extract a subset of columns
                        'af_UKBB', 'ld', 'rsid', 'pos_hg17', 'pos_hg18', 'pos_hg38', '_NUM_ID_')]
colnames(df_beta1)[ match( c('_NUM_ID_'), colnames(df_beta1)) ] = c('_NUM_ID_map_ldref') # change name to avoid warning from combining tables; _NUM_ID_map_ldref is the indices in map_ldref
df_beta1[1:10, ]
dim(df_beta1) 


### Match df_beta1 and map
df_beta2 = snp_match(df_beta1, map, return_flip_and_rev = T)

df_beta2[1:10, ]
dim(df_beta2) 

### Remove missing in G
extracted_sub_G = G[1:nrow(G), df_beta2$'_NUM_ID_'] # extracted subset of G
dim(extracted_sub_G) 

NA_vec = colSums(is.na(extracted_sub_G)) # number of NAs in 950 samples for each snp
which(NA_vec > 0) # No snps have missing values

### QC steps (Reference: https://github.com/privefl/paper-ldpred2/blob/master/code/example-with-provided-ldref.R)
# NOTE: af_UKBB is allele frequency from UKBB
sd_ldref = with(df_beta2, sqrt(2 * af_UKBB * (1 - af_UKBB))) # standard deviations of genotypes of individuals
sd_ss = with(df_beta2, 2 / sqrt(n_eff * beta_se^2)) # standard deviations derived from the summary statistics 
# Note: the numerator was 1 in their online code, but they changed it to 2 not long ago. If using 2, we would remove all snps, so we use 1 here.

is_bad = sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05 
sum(is_bad) 

### Remove bad snps
df_beta_final = df_beta2[!is_bad, ]
extracted_sub_G = extracted_sub_G[, !is_bad]

### Reverse the beta in df_beta_final with flipped alleles so that its value corresponds to map_ldref or reversed G rather than G
reverse_ids = which(df_beta_final$'_REV_') # rows for reverse
df_beta_final$beta[reverse_ids] = - df_beta_final$beta[reverse_ids] # reverse beta
extracted_sub_G[, reverse_ids] = 2 - extracted_sub_G[, reverse_ids] # modify genotype
rownames(extracted_sub_G) = obj_bigSNP$fam$sample.ID # add sample names

dim(df_beta_final) 

dim(extracted_sub_G) 

df_beta_final[1:5, ] # _NUM_ID_ is the ids in extracted_sub_G; _NUM_ID_map_ldref is the ids in map_ldref
extracted_sub_G[1:5, 1:10]


### 3. Find correlations ------------------------------------------------------------------------------------------------------------------------------------------------
tmp = tempfile(tmpdir = "cor_mat")

for (chr in 1:22) {
  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta_final' for current chr
  ind.chr = which(df_beta_final$chr == chr)
  ## indices in 'map_ldref' for current chr
  ind.chr2 = df_beta_final$"_NUM_ID_map_ldref"[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 = match(ind.chr2, which(map_ldref$chr == chr))
  
  corr_chr = readRDS(paste0("LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    corr = as_SFBM(corr_chr, tmp)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}
dim(corr) 

### 4. Find scores ------------------------------------------------------------------------------------------------------------------------------------------------
### Estimate of h2 from LD Score regression
(ldsc = with(df_beta_final, snp_ldsc(ld, 
                                     length(ld), 
                                     chi2 = (beta / beta_se)^2,
                                     sample_size = n_eff, blocks = NULL)))
h2_est = ldsc[["h2"]] 

### LDpred2-auto
set.seed(111)
multi_auto = snp_ldpred2_auto(corr, 
                              df_beta_final, 
                              h2_init = h2_est,
                              vec_p_init = seq_log(1e-4, 0.3, length.out = 30),
                              # allow_jump_sign = FALSE, 
                              # shrink_corr = 0.95,
                              ncores = NCORES)
beta_auto_full = sapply(multi_auto, function(auto) auto$beta_est) # 966,584 x 30
dim(beta_auto_full)
(sapply(multi_auto, function(auto) auto$p_est)) # 30 estimated p
(sapply(multi_auto, function(auto) auto$h2_est)) # 30 estimated h2

mean(sapply(multi_auto, function(auto) auto$p_est)) 
mean(sapply(multi_auto, function(auto) auto$h2_est)) 
mean(sapply(multi_auto, function(auto) auto$beta_est)) 


(range = sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep = (range > (0.95 * quantile(range, 0.95))))
beta_auto = rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto = as.vector(extracted_sub_G %*% beta_auto) # final predictions
names(pred_auto) = obj_bigSNP$fam$sample.ID # add names
