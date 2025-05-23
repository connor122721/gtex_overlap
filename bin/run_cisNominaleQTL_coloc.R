# Connor Murray
# Started 12.8.2024; modified 5.1.25 colocalization w/HF GWAS
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(coloc)
library(tidyverse)
library(foreach)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="Processed GWAS summary statistics.")
parser$add_argument("--eqtl", required=TRUE, help="cis-eQTL statistics.")
parser$add_argument("--shortList", required=TRUE, help="Extracts candidate SNPs/eGenes to test for colocalization from first cis eqtl-tensorQTL run.")
parser$add_argument("--chromosome", required=TRUE, help="Current chromosome.")
parser$add_argument("--N_gwas", required=TRUE, help="Samples in GWAS.")
parser$add_argument("--N_eqtl", required=TRUE, help="Samples in eQTL study.")
parser$add_argument("--prefix", required=TRUE, help="Prefix of output files.")

args <- parser$parse_args()

gwas_input <- args$gwas
eqtl_inpt <- args$eqtl
sig_genes <- unlist(fread(args$shortList) %>% 
                      select(phenotype_id))
chromosome <- args$chromosome
N_gwas <- as.numeric(args$N_gwas)
N_eqtl <- as.numeric(args$N_eqtl)
output_pre <- args$prefix

# setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/");gwas_input="output/gwas/processed_jurgens24_gwas_HF_chr7.rds"
# eqtl_inpt="/standard/vol185/cphg_Manichaikul/projects/GTEx/GTEx_v10_TissueSpecific-eQTL/Heart_Left_Ventricle.v10.allpairs.chr7.parquet";shortlist="output/tensorqtl/topchef_chr7_MaxPC70.cis_qtl.txt.gz"
# N_gwas=1665481;N_eqtl=516;chromosome="chr7"; output_pre="Jurgens2024"
# sig_genes=unlist(fread(shortlist) %>% select(phenotype_id))

### Datasets & Setup ###

# By chromosome coloc
coloc_chrom <- function(chr) {
  
  # Start on chromosome: chr="chr7"
  chromi=chr
  print(paste("Running:", chromi, sep=" "))
  
  # Read in parquet files of GTEx
  qtl <- data.table(arrow::read_parquet(eqtl_inpt)) %>% 
    mutate(chrom = chromi,
           snp = as.numeric(tstrsplit(variant_id, "_")[[2]]),
           variant_id = paste(chrom, snp, sep=":"),
           type = "quant",
           maf = case_when(af > 0.5 ~ 1-af,
                           TRUE ~ af),
           n = round(ma_count/(2*maf)),
           phenotype_id = tstrsplit(gene_id, ".", fixed=T)[[1]])
  
  # Read in GWAS
  gwas <- data.table(readRDS(gwas_input))
  
  # Progress
  print(paste("Done reading/preparing data:", chromi))
  
    # Window size coloc function between Shah and eQTL datasets
    colocWindow <- function(dt1, N1=N_eqtl, focalGene, dt2, N2=N_gwas) {
      # dt1=qtl; dt2=gwas; focalGene="ENSG00000128591"; N1=516; N2=1665481
      
      # Restrict to all eQTL SNPs 
      dt1 <- dt1[phenotype_id %in% focalGene][!duplicated(variant_id)]
      min_window <- as.numeric(min(dt1$snp))
      max_window <- as.numeric(max(dt1$snp))
      
      # Ensure same length
      dt2 <- dt2[snpID_hg38 %in% dt1$variant_id]
      dt1 <- dt1[variant_id %in% dt2$snpID_hg38]
      
      # Combine
      dtsub <- dt1 %>% 
        left_join(dt2 %>% 
              select(variant_id=snpID_hg38,
                     beta_gwas=beta, 
                     varbeta_gwas=standard_error,
                     pvalue_gwas=p_value,
                     A1,
                     A2,
                     rsid=SNP))
      
      # Number of putative significant cis regions
      nSnps = as.numeric(dtsub %>% 
              filter(pval_nominal <= 1e-5) %>% 
              summarize(n()))
      
      # Make list for colocalization (TOPChef)
      dt1 <- list(snp = dt1$variant_id,
                  position = dt1$snp,
                  type = "quant",
                  N = N1,
                  MAF = as.numeric(dt1$maf),
                  pvalues = as.numeric(dt1$pval_nominal),
                  beta = as.numeric(dt1$slope),
                  varbeta = as.numeric(dt1$slope_se)^2)
      
      # Make list for colocalization (Shah)
      dt2 <- list(snp = dt2$snpID_hg38,
                  position = dt2$pos_hg38,
                  type = "quant",
                  N = N2,
                  MAF = as.numeric(dt2$maf),
                  pvalues = as.numeric(dt2$p_value),
                  beta = as.numeric(dt2$beta),
                  varbeta = as.numeric(dt2$standard_error)^2)
      
      # Colocalization analysis using the coloc package
      coloc_res <- coloc.abf(dataset1 = dt1,
                             dataset2 = dt2)
      
      maxi = max(coloc_res$results$SNP.PP.H4)
      snpi = coloc_res$results[which(coloc_res$results$SNP.PP.H4==maxi),]$snp
      
      # Find candidate sentinenal variant
      dtsub[, rank_eqtl := rank(pval_nominal), by = phenotype_id]
      dtsub[, rank_gwas := rank(pvalue_gwas), by = phenotype_id]
      dtsub[, total_rank := rank_eqtl + rank_gwas]
      
      # For each gene, choose the variant with the smallest total_rank
      sentinel_combined <- dtsub[order(rank_gwas), .SD[1], by = phenotype_id]
      sent = sentinel_combined$variant_id
      
      #plot(coloc_res)
      
      # Results of colocalization
      co <- data.table(chrom=chromi,
                       minPos=min_window,
                       maxPos=max_window,
                       gene=focalGene,
                       min_p.eqtl=min(dt1$pvalues),
                       min_p.gwas=min(dt2$pvalues),
                       nsnps=coloc_res$summary["nsnps"],
                       nCis=nSnps,
                       PP.H0=coloc_res$summary["PP.H0.abf"],
                       PP.H1=coloc_res$summary["PP.H1.abf"],
                       PP.H2=coloc_res$summary["PP.H2.abf"],
                       PP.H3=coloc_res$summary["PP.H3.abf"],
                       PP.H4=coloc_res$summary["PP.H4.abf"],
                       H4_H3_ratio=coloc_res$summary["PP.H4.abf"]/coloc_res$summary["PP.H3.abf"],
                       sentinel=sent,
                       sent_maf=sentinel_combined$maf,
                       sent_dist_TSS=sentinel_combined$start_distance,
                       beta_qtl=sentinel_combined$slope,
                       betase_qtl=sentinel_combined$slope_se,
                       beta_gwas=sentinel_combined$beta_gwas,
                       betase_gwas=sentinel_combined$varbeta_gwas,
                       A1=sentinel_combined$A1,
                       A2=sentinel_combined$A2,
                       maxSNP=snpi,
                       maxPP.H4=maxi,
                       gwas_pre=output_pre)
      
      # Finish
      return(co)
    }
    
  # Run coloc analysis with error handling
  coloc_results <- map(sig_genes, possibly(function(gene) {
    colocWindow(dt1 = qtl, 
                dt2 = gwas, 
                focalGene = gene)}, otherwise = NULL))
  
  # Filter out NULL results and keep only data frames
  coloc_results_valid <- keep(coloc_results, is.data.frame)
  
  # Combine results with rbindlist
  fin_coloc <- rbindlist(coloc_results_valid, fill = TRUE)
  
  # Finish
  print(paste("Finish:", chr))
  return(fin_coloc)

}

# Run code
dt <- coloc_chrom(chromosome)

# Output results
write_delim(dt, file = paste("coloc_eqtl_", output_pre, "_", chromosome, ".txt", sep=""), delim = "\t")
