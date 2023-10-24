# EDF6: Forest plots of genetic associations for selected traits across different ancestry groups
# beta_estimate_sysmex.xlsx contains beta estimates from the top GWAS hits under perturbation conditions

library(ggplot2)
library(glue)
library(data.table)
library(parallel)
library(metafor)
library(dplyr)

sysmex_estimates <- data.table(readxl::read_excel('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/beta_estimate_sysmex.xlsx'))
sysmex_estimates$trait <- gsub('-q', '', sysmex_estimates$trait)

sysmex_estimates <- sysmex_estimates[traitname %in% c(
  "ret|Rotenone 6h overnight|RET1_CV_SFL-q", "wdf|KCl 17h|NE2_NE4_ratio-q",
  "wdf|LPS 18h|NE4_SD_SFL-q", "wdf|Alhydrogel 21h|NE4_SD_SFL-q",
  "wdf|Colchicine 20h|NE4_SD_SFL-q", "wdf|Pam3CSK4 19h|NE1_Med_FSC-q" 
                              )]

all_results <-
mclapply(c('asian', 'black', 'hispanic', 'other'), function(race){
  
  sysmex_estimates_anc <- sysmex_estimates
  sysmex_estimates_anc$gwas_paths <- glue('"/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/gwas/{sysmex_estimates$channel}/perturbations/{sysmex_estimates$ptb_name}/non_eur/maf001_{race}.{sysmex_estimates$trait}-q.glm.linear.zst"')
  
  results <-
    lapply(1:nrow(sysmex_estimates_anc), function(i){
      input_path <- sysmex_estimates_anc[i,]$gwas_paths
      sumstats <- fread(cmd=paste('zstdcat', input_path))
      
      # match using chrom and pos
      res <- sumstats[(`#CHROM` == sysmex_estimates_anc[i,]$CHR) & (POS == sysmex_estimates_anc[i,]$POS),]
      res$gene <- sysmex_estimates_anc[i,]$gene
      res$trait <- sysmex_estimates_anc[i,]$trait
      res$traitname <- as.character(glue('{sysmex_estimates_anc[i,]$channel}|{sysmex_estimates_anc[i,]$ptb_name}|{sysmex_estimates_anc[i,]$trait}'))
      res$ancestry <- race
      res
    })
  results <- rbindlist(results)
  results
})

all_results <- rbindlist(all_results)
sysmex_estimates$ancestry <- 'european'

all_sysmex_results <- rbind(
  sysmex_estimates[, c('gene', 'BETA', 'SE', 'P', 'trait', 'traitname', 'ancestry', 'ID', 'REF', 'ALT', 'A1', 'OBS_CT')],
  all_results[, c('gene', 'BETA', 'SE', 'P', 'trait', 'traitname', 'ancestry', 'ID', 'REF', 'ALT', 'A1', 'OBS_CT')])

all_sysmex_results$abs_beta <- abs(all_sysmex_results$BETA)

# if A1 is not the ALT allele, then it is flipped relative to the EUR ancestry
all_sysmex_results$beta_aligned <- all_sysmex_results$BETA
all_sysmex_results[A1 != ALT]$beta_aligned <- -1*all_sysmex_results[A1 != ALT]$BETA

ggplot(all_sysmex_results, aes(x=beta_aligned, y=gene, shape=ancestry, color=ancestry)) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin =beta_aligned - SE, xmax =beta_aligned + SE), height = 0.2) +
  labs(x = "Beta", y = "Gene") +
  theme_minimal() +
  coord_flip()
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/cross_ancestry.pdf', width = 6, height = 6)
fwrite(all_sysmex_results, '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/beta_estimates_cross_ancestry.tsv', sep='\t')

# recode ancestry names
ancestry_names <- c('european' = 'EUR', 'black' = 'AFR', 'hispanic' = 'HISP', 'asian'='ASIAN', 'other'='OTHER')
all_sysmex_results$ancestry <- recode(all_sysmex_results$ancestry, !!!ancestry_names)
all_sysmex_results[['ancestry_count']] <- paste0(all_sysmex_results$ancestry, ' (', all_sysmex_results$OBS_CT, ')')
all_sysmex_results <- all_sysmex_results[ancestry %in% c('EUR', 'ASIAN', 'AFR', 'OTHER')]
all_sysmex_results <- all_sysmex_results[gene %in% c("RHCE", "BCL2A1", "HK1", "TLR1")]

unique_genes <- unique(all_sysmex_results$gene)
results_list <- list()
plot_list <- list()

pdf('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/cross_ancestry_metafor.pdf', width = 8, height = 6)
par(mfrow=c(2, 2))
for (gene_ in unique_genes) {
  subset_data <- all_sysmex_results[gene == gene_]
  res <- rma(yi=beta_aligned, sei=SE, data=subset_data, method="REML")
  results_list[[gene_]] <- res
  title <- paste0(gene_, '\n', gsub('\\|', ' ', gsub('_', ' ', gsub('-q', '', subset_data$traitname[[1]]))))
  forest(res, slab=subset_data$ancestry_count, main=title)
}
dev.off()
