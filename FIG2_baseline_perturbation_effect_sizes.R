# Figure 2B
# Compare the beta coefficients of select traits under perturbation with baseline and with prior GWAS studies (from GWAS Catalog)
# Figure 2A is based on a customized version of Fujiplot: https://github.com/mkanai/fujiplot

library(ggplot2)
library(glue)

sysmex_estimates <- data.table(readxl::read_excel('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/beta_estimate_sysmex.xlsx'))
sysmex_estimates$trait <- gsub('-q', '', sysmex_estimates$trait)
sysmex_estimates_baseline <- sysmex_estimates
sysmex_estimates_baseline$baseline_paths <- glue('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/gwas/{sysmex_estimates$channel}/perturbations/Baseline/maf005_white.{sysmex_estimates$trait}-q.glm.linear.zst.ADD.gz')

baseline_results <-
lapply(1:nrow(sysmex_estimates_baseline), function(i){
  sumstats <- fread(sysmex_estimates_baseline[i,]$baseline_paths)
  res <- sumstats[ID == sysmex_estimates_baseline[i,]$ID]
  res$gene <- sysmex_estimates_baseline[i,]$gene
  res$trait <- sysmex_estimates_baseline[i,]$trait
  res$traitname <- glue('{sysmex_estimates_baseline[i,]$channel}|Baseline|{sysmex_estimates_baseline[i,]$trait}')
  res
})
baseline_results <- rbindlist(baseline_results)
baseline_results$traitname <- as.character(baseline_results$traitname)

sysmex_estimates$type <- 'Perturbation'
baseline_results$type <- 'Baseline'

all_sysmex_results <- rbind(sysmex_estimates[, c('gene', 'BETA', 'SE', 'P', 'trait', 'traitname', 'type')],
                            baseline_results[, c('gene', 'BETA', 'SE', 'P', 'trait', 'traitname', 'type')])

chen_estimates <- data.table(readxl::read_excel('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/beta_estimates_blood_traits.xlsx'))
chen_estimates <- chen_estimates[, c('gene', 'beta_abs', 'se', 'pValue', 'traitName', 'author')]
chen_estimates$trait <- chen_estimates$traitName
setnames(chen_estimates, 'author', 'type')

all_sysmex_results$BETA <- abs(all_sysmex_results$BETA)
setnames(all_sysmex_results, c('gene', 'beta_abs', 'se', 'pValue', 'trait', 'traitName', 'type'))

beta_estimates <-
rbind(
  chen_estimates,
  all_sysmex_results
)

# figure was generated using Prism using the table below
fwrite(beta_estimates, '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/beta_estimates_joint.tsv', sep='\t')

# additional plots (not used)
ggplot(beta_estimates, aes(x=beta_abs, y=gene, shape=type, color=trait)) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin = beta_abs - se, xmax = beta_abs + se), height = 0.2) +
  labs(x = "Absolute Beta", y = "Gene") +
  theme_minimal() +
  coord_flip()

ggplot(beta_estimates, aes(x=beta_abs, y=trait, shape=type, color=trait)) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin = beta_abs - se, xmax = beta_abs + se), height = 0.8) +
  labs(x = "Absolute Beta", y = "Gene") +
  theme_minimal() +
  scale_y_discrete(breaks = NULL) +
  theme(legend.position="bottom") +
  facet_wrap('.~gene', scales='free_x') +
  coord_flip()
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/beta_comparisons/absolute_beta.pdf', width = 7, height = 7)  
