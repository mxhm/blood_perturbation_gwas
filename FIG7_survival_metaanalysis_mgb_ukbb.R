# Figure 7A, 7B, 7C: plots for MGB and UKBB PRS-disease associations

# this script reads in the results from MGB and UKBB, does a meta-analysis
# plots a subset of the associations for the main figure 7B
# calculates ICA projection of the PRS-disease assocations in figure 7C
# the panels in 7A are based on
# - revision/survival_subset_ukbb_230706
# - revision/survival_subset_mgb_230706
# - metanalyses: revision/plots/meta_analyses_counting_230706

library(data.table)
library(metafor)
library(qvalue)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(parallel)
library(stringr)
library(colorspace)
library(dendsort)
library(fastICA)
library(ggrepel)

# meta_stats_230706.tsv is renamed to S6_PGS_meta-analysis.tsv
mgb_ukbb_traits <- data.table(readxl::read_excel('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/clinical_traits_icd10_mgb_ukbb_matching.xlsx'))

mgb_survival_pvals_sex_pcs <-
  rbindlist(lapply(
    list.files(
      path="/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/survival_ct_p01_sex_pcs2_cox_stratsex_mgb230703/",
      pattern="*.stats.txt",
      full.names=TRUE,
    ),
    fread))


ukbb_survival_pvals_sex_pcs <-
  rbindlist(lapply(
    list.files(
      path="/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/survival_ct_p01_sex_pc2_counting_ukbb230703/",
      pattern="*.stats.txt",
      full.names=TRUE,
    ),
    fread))


#### debug plots: comparing fit3 and fit5 outputs (delayed entry vs entry at birth) ####
plot_cols <- c('clinical_trait', 'blood_trait', 'term', 'statistic', 'estimate', 'method', 'p.value')

wide_dt <- dcast(mgb_survival_pvals_sex_pcs[method %in% c('fit5', 'fit4', 'fit3') & (term == 'risk_score'), ..plot_cols],
                 clinical_trait + blood_trait + term ~ method, value.var = c('statistic', 'estimate', 'p.value'))

# calculate the -log10 values
wide_dt$log_pvalue_fit3 <- -log10(wide_dt$p.value_fit3)
wide_dt$log_pvalue_fit5 <- -log10(wide_dt$p.value_fit5)

# fit a linear model
fit <- lm(log_pvalue_fit5 ~ log_pvalue_fit3, data = wide_dt)

# get the slope and intercept
slope <- coef(fit)[2]
intercept <- coef(fit)[1]

# calculate the correlation coefficient and p-value
cor_test <- cor.test(wide_dt$log_pvalue_fit3, wide_dt$log_pvalue_fit5)

cor_coeff <- cor_test$estimate
p_value <- cor_test$p.value

list(slope = slope, intercept = intercept, correlation_coefficient = cor_coeff, p_value = p_value)
     

ggplot(wide_dt, aes(-log(p.value_fit3), -log(p.value_fit5))) +
  geom_point(aes(size=-log(p.value_fit3)), alpha=0.1) +
  geom_smooth(method = 'lm', se = TRUE, color = 'blue', fill = 'lightblue') +
  geom_abline(slope = 1, intercept = 0, color='red') +
  theme(legend.position="bottom")

ggplot(wide_dt, aes(estimate_fit3, estimate_fit5)) +
  geom_point(aes(size=-log(p.value_fit3)), alpha=0.1) +
  geom_smooth(method = 'lm', se = TRUE, color = 'blue', fill = 'lightblue') +
  geom_abline(slope = 1, intercept = 0, color='red') +
  theme(legend.position="bottom")

ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/lateentry_vs_normal_cox_estimates.pdf', width = 8, height = 9)

#### data preprocessing
ukbb_survival_pvals_sex_pcs <- ukbb_survival_pvals_sex_pcs[method == 'fit3' & term == 'risk_score']
mgb_survival_pvals_sex_pcs <- mgb_survival_pvals_sex_pcs[method == 'fit3' & term == 'risk_score']

# p-value adjustment is done using the meta-analysis output
ukbb_sel <- ukbb_survival_pvals_sex_pcs
mgb_sel <- mgb_survival_pvals_sex_pcs

mgb_ukbb_traits$UKBB2 <- gsub('date_', 'age_', mgb_ukbb_traits$UKBB)
mgb_sel$name <- mgb_sel$clinical_trait

ukbb_sel2 <- merge(ukbb_sel, mgb_ukbb_traits, by.x = 'clinical_trait', by.y = 'UKBB2')
unique(mgb_sel$name)
unique(ukbb_sel2$name)
unique(mgb_ukbb_traits$UKBB)
unique(mgb_ukbb_traits$UKBB2)

joint_assoc <- merge(ukbb_sel2, mgb_sel, by= c('blood_trait', 'name'))
mgb_sel$site <- 'MGB'
ukbb_sel2$site <- 'UKBB'

meta_stats_all <-
mclapply(1:nrow(joint_assoc), mc.cores = 16, function(i){
    
  name_ <- joint_assoc[i,]$name
  blood_trait_ <- joint_assoc[i,]$blood_trait
  
  clin_trait_ukbb <- joint_assoc[i,]$clinical_trait.x
  clin_trait_mgb <- joint_assoc[i,]$clinical_trait.y
  
  joint_stats <-
    rbind(
      ukbb_sel2[(blood_trait == blood_trait_) & (clinical_trait == clin_trait_ukbb)],
      mgb_sel[(blood_trait == blood_trait_) & (clinical_trait == clin_trait_mgb)],
      fill=TRUE
    )
  
  if (nrow(joint_stats) > 2){
    print('too many rows:')
    print(joint_stats)
  }
  
  # Perform the meta-analysis
  meta_res <- rma(yi = estimate, sei = std.error, data = joint_stats, slab=paste0(site, ' ', substr(clinical_trait, 1, 15)))
  meta_stats_df <- data.frame(broom::tidy(meta_res))
  meta_stats_df$blood_trait <- blood_trait_
  meta_stats_df$name <- name_
  meta_stats_df$ukbb <- clin_trait_ukbb
  meta_stats_df$mgb <- clin_trait_mgb
  meta_stats_df$mgb_raw <- joint_stats[site == 'MGB',]$p.value
  meta_stats_df$ukbb_raw <- joint_stats[site == 'UKBB',]$p.value
  
  if (meta_stats_df$p.value < 0.01 & meta_stats_df$mgb_raw < 0.05 & meta_stats_df$ukbb_raw < 0.05){
    # Print the summary
    print(paste0(name_, ' - ', blood_trait_))
    print(summary(meta_res))
    pdf(paste0('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/meta_analyses_counting_230706/', blood_trait_, '_', clin_trait_mgb, '-', clin_trait_ukbb, '.pdf'))
    forest(meta_res,
           header = paste0(blood_trait_, '\n', name_),
           )
    dev.off()
  }
  meta_stats_df
}
)

meta_stats_all <- rbindlist(meta_stats_all)
hist(meta_stats_all$p.value)

meta_stats_all$meta_qval <- qvalue(meta_stats_all$p.value)$qvalue
meta_stats_all$mgb_qval <- qvalue(meta_stats_all$mgb_raw)$qvalue
meta_stats_all$ukbb_qval <- qvalue(meta_stats_all$ukbb_raw)$qvalue

fwrite(meta_stats_all, '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/meta_stats_230706.tsv', sep='\t')
fwrite(meta_stats_all[meta_qval < 0.1,], '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/meta_stats_q0.1_230706.tsv', sep='\t')

# plot a subset in main figure 7B, full dataset in supplement
meta_stats_all_plot <- meta_stats_all[(meta_qval < 0.05) & (mgb_qval < 0.25) & (ukbb_qval < 0.2),]
unique(meta_stats_all_plot$name)

length(unique(meta_stats_all$name))
length(unique(meta_stats_all$blood_trait))

#### heatmaps of meta-analysis results

# make sure none of the UKBB-MGB pairs are duplicates
meta_stats_all_dedup <- meta_stats_all

beta_coeff <- data.frame(dcast(meta_stats_all_dedup, blood_trait ~ name, value.var = 'estimate'))
qvalues <- data.frame(dcast(meta_stats_all_dedup, blood_trait ~ name, value.var = 'meta_qval'))
tscores <- data.frame(dcast(meta_stats_all_dedup, blood_trait ~ name, value.var = 'statistic'))
rownames(beta_coeff) <- beta_coeff$blood_trait
rownames(qvalues) <- beta_coeff$blood_trait
rownames(tscores) <- beta_coeff$blood_trait

sum(is.na(tscores))

tscores <- tscores[,-1]
tscores[is.na(tscores)] <- 0
beta_coeff[is.na(beta_coeff)] <- 0


# plot only associations with qvalue <= 0.1
beta_coeff_zeroed <- beta_coeff[, -1 ]
beta_coeff_zeroed[(qvalues[, -1] >= 0.05) |
                        (is.na(qvalues[, -1]))] <- 0

qvals <- qvalues[, -1]
qvals[is.na(qvals)] <- 1

qval_labels <- matrix(nrow=nrow(qvals), ncol=ncol(qvals))
qval_labels[qvals <= 1] <- ''
qval_labels[qvals <= 0.05] <- '·'
qval_labels[qvals <= 0.01] <- '··'
qval_labels[qvals <= 0.001] <- '···'

colnames(qval_labels) <- colnames(qvalues)[-1]
rownames(qval_labels) <- rownames(qvalues)

sign_blood_traits <- rownames(beta_coeff)
sign_clinical_traits <- colnames(beta_coeff[,-1])

mgb_ukbb_traits$name2 <- make.names(mgb_ukbb_traits$name)

clinical_annotation_df <- data.frame(na.exclude(mgb_ukbb_traits[, c('name2', 'category')]))
row.names(clinical_annotation_df) <- clinical_annotation_df$name2

# TODO create annotations from the readable blood trait names
bloodtrait_annotation_df <- meta_stats_all_dedup[, c('blood_trait')]
bloodtrait_annotation_df <- as.data.frame(unique(bloodtrait_annotation_df))
rownames(bloodtrait_annotation_df) <- bloodtrait_annotation_df$blood_trait

bloodtrait_annotation_df$blood_trait_original <- bloodtrait_annotation_df$blood_trait
bloodtrait_annotation_df$blood_trait_original <- 
  gsub('Baseline', 'Baseline_0h', bloodtrait_annotation_df$blood_trait_original)
bloodtrait_annotation_df$blood_trait_original <- 
  gsub('-3', '3', bloodtrait_annotation_df$blood_trait_original)
bloodtrait_annotation_df$blood_trait_original <- 
  gsub('overnight_', '', bloodtrait_annotation_df$blood_trait_original)

bloodtrait_annotation_df$celltype <- str_split_fixed(bloodtrait_annotation_df$blood_trait_original, '-', n=3)[,3]
bloodtrait_annotation_df$celltype <- str_split_fixed(bloodtrait_annotation_df$celltype, '_', n=2)[,1]
bloodtrait_annotation_df$celltype <- gsub('[0-9]', '', bloodtrait_annotation_df$celltype)
bloodtrait_annotation_df$channel <- str_split_fixed(bloodtrait_annotation_df$blood_trait_original, '-', n=2)[,1]

# define colors for the annotation
channel_cols <- brewer.pal(n = 4, name = "Set1")
names(channel_cols) <- unique(bloodtrait_annotation_df$channel)
celltype_cols <- darken(brewer.pal(n = length(unique(bloodtrait_annotation_df$celltype)), name = "Set3"), 0.05)
names(celltype_cols) <- unique(bloodtrait_annotation_df$celltype)
category_cols <- darken(brewer.pal(n = length(unique(clinical_annotation_df$category)), name = "Spectral"), 0.05)
names(category_cols) <- unique(clinical_annotation_df$category)

category_cols_dark <- darken(brewer.pal(n = length(unique(clinical_annotation_df$category)), name = "Spectral"), 0.25)
names(category_cols_dark) <- unique(clinical_annotation_df$category)

ann_colors = list(
  channel = channel_cols,
  celltype = celltype_cols,
  category = category_cols,
  category_cols_dark = category_cols_dark
)
#saveRDS(ann_colors, file='/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/ann_colors_prs.rds')


 
callback = function(hc, ...){dendsort(hc)}

rownames(beta_coeff) <- gsub('_', ' ', rownames(beta_coeff))
rownames(beta_coeff) <- gsub('-', ' ', rownames(beta_coeff))

rownames(qval_labels) <- gsub('_', ' ', rownames(qval_labels))
rownames(qval_labels) <- gsub('-', ' ', rownames(qval_labels))

rownames(bloodtrait_annotation_df) <- gsub('_', ' ', rownames(bloodtrait_annotation_df))
rownames(bloodtrait_annotation_df) <- gsub('-', ' ', rownames(bloodtrait_annotation_df))

colnames(beta_coeff) <- gsub('\\.', ' ', tolower(colnames(beta_coeff)))
colnames(qval_labels) <- gsub('\\.', ' ', tolower(colnames(qval_labels)))
rownames(clinical_annotation_df) <- gsub('\\.', ' ', tolower(rownames(clinical_annotation_df)))

rg <- 0.3
pheatmap(
  beta_coeff[,-1],
  display_numbers = qval_labels,
  color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
  breaks = seq(-rg, rg, length.out = 100),
  cutree_rows=7,
  cutree_cols=7,
  #clustering_distance_rows='correlation',
  clustering_distance_cols='correlation',
  annotation_col=clinical_annotation_df[,2, drop=FALSE],
  annotation_row = bloodtrait_annotation_df[, c('channel', 'celltype'), drop=FALSE],
  annotation_colors = ann_colors,
  fontsize_row = 3,
  fontsize_number = 14,
  border_color = NA,
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/plink_ct_p01_meta_all_230706.pdf',
  width = 12,
  height = 20,
  angle_col="45",
  treeheight_row=30,
  treeheight_col=30,
  clustering_callback = callback
)

#### make a plot of significant association subset ####
sign_blood <- gsub('[_-]', ' ', unique(meta_stats_all_dedup[meta_qval < 0.05]$blood_trait))
sign_clin <- tolower(gsub('\\.', ' ', make.names(unique(meta_stats_all_dedup[meta_qval < 0.05]$name))))

rg <- 0.3
pheatmap(
  beta_coeff[sign_blood, sign_clin],
  display_numbers = qval_labels[sign_blood, sign_clin],
  color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
  breaks = seq(-rg, rg, length.out = 100),
  cutree_rows=6,
  cutree_cols=6,
  #clustering_distance_rows='correlation',
  clustering_distance_cols='correlation',
  annotation_col=clinical_annotation_df[,2, drop=FALSE],
  annotation_row = bloodtrait_annotation_df[, c('channel', 'celltype'), drop=FALSE],
  annotation_colors = ann_colors,
  fontsize_row = 6,
  fontsize_number = 8,
  border_color = NA,
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/plink_ct_p01_meta_subset_230706.pdf',
  width = 9,
  height = 14,
  angle_col="45",
  treeheight_row=30,
  treeheight_col=30
)

#### plot subset where the meta-analysis is significant, and the two cohorts qval <0.25
meta_stats_all_dedup$meta_qval <- qvalue(meta_stats_all_dedup$p.value)$qvalue
meta_stats_all_dedup$mgb_qval <- qvalue(meta_stats_all_dedup$mgb_raw)$qvalue
meta_stats_all_dedup$ukbb_qval <- qvalue(meta_stats_all_dedup$ukbb_raw)$qvalue

meta_stats_all_dedup_plot <- meta_stats_all_dedup[(meta_qval < 0.05) & (mgb_qval < 0.25) & (ukbb_qval < 0.2),]

sign_blood <- meta_stats_all_dedup_plot$blood_trait
sign_clin <- meta_stats_all_dedup_plot$name

meta_stats_all_subset <- meta_stats_all_dedup[
  (blood_trait %in% sign_blood) & (name %in% sign_clin)]

beta_coeff <- data.frame(dcast(meta_stats_all_subset, blood_trait ~ name, value.var = 'estimate'))
qvalues <- data.frame(dcast(meta_stats_all_subset, blood_trait ~ name, value.var = 'meta_qval'))
tscores <- data.frame(dcast(meta_stats_all_subset, blood_trait ~ name, value.var = 'statistic'))
rownames(beta_coeff) <- beta_coeff$blood_trait
rownames(qvalues) <- beta_coeff$blood_trait
rownames(tscores) <- beta_coeff$blood_trait

tscores <- tscores[,-1]
tscores[is.na(tscores)] <- 0
beta_coeff[is.na(beta_coeff)] <- 0

qvals <- qvalues[, -1]
qvals[is.na(qvals)] <- 1

qval_labels <- matrix(nrow=nrow(qvals), ncol=ncol(qvals))
qval_labels[qvals <= 1] <- ''
qval_labels[qvals <= 0.05] <- '·'
qval_labels[qvals <= 0.01] <- '··'
qval_labels[qvals <= 0.001] <- '···'

colnames(qval_labels) <- colnames(qvalues)[-1]
rownames(qval_labels) <- rownames(qvalues)

rownames(beta_coeff) <- gsub('_', ' ', rownames(beta_coeff))
rownames(beta_coeff) <- gsub('-', ' ', rownames(beta_coeff))

rownames(qval_labels) <- gsub('_', ' ', rownames(qval_labels))
rownames(qval_labels) <- gsub('-', ' ', rownames(qval_labels))

rownames(bloodtrait_annotation_df) <- gsub('_', ' ', rownames(bloodtrait_annotation_df))
rownames(bloodtrait_annotation_df) <- gsub('-', ' ', rownames(bloodtrait_annotation_df))

rg <- 0.1
pheatmap(
  beta_coeff[,-1],
  display_numbers = qval_labels,
  color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
  breaks = seq(-rg, rg, length.out = 100),
  cutree_rows=5,
  cutree_cols=3,
  #clustering_distance_rows='correlation',
  clustering_distance_cols='correlation',
  annotation_col=clinical_annotation_df[,2, drop=FALSE],
  annotation_row = bloodtrait_annotation_df[, c('channel', 'celltype'), drop=FALSE],
  annotation_colors = ann_colors,
  fontsize_row = 8,
  fontsize_number = 8,
  border_color = NA,
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/plink_ct_p01_meta_subset2_230706.pdf',
  width = 6.5,
  height = 8,
  angle_col="45",
  treeheight_row=20,
  treeheight_col=20,
  clustering_callback = callback
)

### Figure 7C: ICA projection
sign_blood <- meta_stats_all_dedup[meta_qval < 0.05]$blood_trait
sign_clin <- meta_stats_all_dedup[meta_qval < 0.05]$name

# project all diseases with significant associations, but using all traits
meta_stats_all_subset <- meta_stats_all_dedup[
   (name %in% sign_clin)]

beta_coeff <- data.frame(dcast(meta_stats_all_subset, blood_trait ~ name, value.var = 'estimate'))
qvalues <- data.frame(dcast(meta_stats_all_subset, blood_trait ~ name, value.var = 'meta_qval'))
tscores <- data.frame(dcast(meta_stats_all_subset, blood_trait ~ name, value.var = 'statistic'))
rownames(beta_coeff) <- beta_coeff$blood_trait
rownames(qvalues) <- beta_coeff$blood_trait
rownames(tscores) <- beta_coeff$blood_trait

tscores <- tscores[,-1]
tscores[is.na(tscores)] <- 0
beta_coeff[is.na(beta_coeff)] <- 0

set.seed(42)
mat <- t(tscores)
ica_result <- fastICA(mat, n.comp = 2, alg.typ = 'parallel')

# get the mixing matrix
S <- ica_result$S

# create a data frame with the first two independent components
df <- data.frame(IC1 = S[, 1], IC2 = S[, 2], Disease = rownames(mat))

# get the mixing matrix
mixing_matrix <- t(ica_result$A)
rownames(mixing_matrix) <- colnames(mat)

# use the traits that were used for detail plots as loadings
sel_loadings <- c(
  "ret-LPS_18h-RBC2_Med_SSC",
  "wnr-Water_15h-WBC2_Med_FSC",
  "wnr-KCl_17h-BASO_Med_SSC",
  "pltf-LPS_18h-PLT_F_Med_SFL",
  "pltf-Pam3CSK4_19h-PLT_F_Percentage",
  "wdf-Alhydrogel_21h-NE2_NE4_ratio"
                  )

# create a data frame with the top traits
df_traits <- data.frame(IC1 = mixing_matrix[, 1], IC2 = mixing_matrix[, 2], Trait = rownames(mixing_matrix))
df_traits <- df_traits[df_traits$Trait %in% sel_loadings,]

# plot the independent components
p <- ggplot(df, aes(x = IC1, y = IC2)) +
  geom_point() +
  geom_text_repel(aes(label = Disease), segment.colour = NA) +
  labs(x = "IC1", y = "IC2") +
  geom_text_repel(data = df_traits, aes(x = IC1, y = IC2, label = Trait), color = 'red') +
  geom_segment(data = df_traits, 
               aes(x = 0, y = 0, xend = IC1, yend = IC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 'red', 
               alpha = 0.5) +
  theme_minimal()

p
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/plink_ct_p01_meta_230706_ica_black.pdf',
       width = 10, height = 10)

df <- merge(df, clinical_annotation_df, by.x = 'Disease', by.y = 'name2')
p <- ggplot(df, aes(x = IC1, y = IC2)) +
  geom_point(aes(color=category)) +
  geom_text_repel(aes(label = Disease, color=category), segment.colour = NA) +
  labs(x = "IC1", y = "IC2") +
  geom_text_repel(data = df_traits, aes(x = IC1, y = IC2, label = Trait), color = 'black') +
  geom_segment(data = df_traits, 
               aes(x = 0, y = 0, xend = IC1, yend = IC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 'black', 
               alpha = 0.5) +
  scale_color_manual(values=category_cols_dark) +
  theme_minimal()

p
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/plink_ct_p01_meta_230706_ica_colors.pdf',
       width = 10, height = 10)


#### plot the survival curves for the selected subset of blood traits above ####

fwrite(meta_stats_all, '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/meta_stats_all_230706.tsv', sep='\t')
fwrite(meta_stats_all[
  (blood_trait %in% sel_loadings & meta_qval < 0.05),], '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/meta_stats_plot_230706.tsv', sep='\t')
