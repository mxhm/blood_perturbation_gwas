# Figure 3: associations of raw blood readouts with clinical endpoints
# plots a heatmap of measured blood trait vs. icd10 code + lab value associations
# creates an ICA projection of blood traits and clinical outcomes
# plots distribution of raw readouts against outcomes
#

library(ggplot2)
library(patchwork)
library(RColorBrewer)

library(corrplot)
library(data.table)
library(Cairo)
library(qvalue)
library(stringr)
library(pheatmap)
library(dendsort)
library(readxl)
library(ggrepel)
library(fastICA)
library(dendsort)
library(ggdist)
library(ggbeeswarm)

#### input data preparation ####
perturbation_details <- 
  read_excel('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/perturbation.time.dose.ID.summary.xlsx')

mgb_ukbb_traits <- data.table(readxl::read_excel('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/clinical_traits_icd10_mgb_ukbb_matching.xlsx'))
lab_trait_annotations <- fread('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/lab_trait_annotations.tsv')
lab_trait_annotations[is.na(clinical_trait)]$clinical_trait <- 'NA'

mgb_ukbb_traits$name <- tolower(mgb_ukbb_traits$name)
mgb_ukbb_traits$name2 <- make.names(tolower(mgb_ukbb_traits$name))
mgb_ukbb_traits$name <- gsub('\\.', ' ', mgb_ukbb_traits$name2)

clinical_annotation_df <- data.frame(na.exclude(mgb_ukbb_traits[, c('name', 'name2', 'category')]))
row.names(clinical_annotation_df) <- clinical_annotation_df$name

clinical_annotation_df$type <- 'diagnosis'
colnames(clinical_annotation_df)
setnames(lab_trait_annotations, c('type', 'name2', 'name', 'category'))
col_order <- c("name","name2","category","type")
lab_trait_annotations <- lab_trait_annotations[, ..col_order]

clinical_annotation_df <-
  rbind(clinical_annotation_df, lab_trait_annotations)
rownames(clinical_annotation_df) <- clinical_annotation_df$name

lab_associations <-
  fread(
    '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/clumps_gated_quantile_lab_regression_notransplant_21-07-06_quantile.csv'
  )
# fix conversion of NA lab to automatic NaN value
lab_associations[is.na(lab)]$lab <- 'NA'

disease_associations <-
  fread(
    '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/clumps_gated_quantile_disease_regression_notransplant_23-07-12_quantile.csv'
  )

lab_associations$type = 'lab'
disease_associations$type = 'disease'

setnames(disease_associations, 'disease', 'clinical_trait')
setnames(lab_associations, 'lab', 'clinical_trait')

# require at least 40 cases in each condition, and drop ret isobutyric readout which has numeric issues
disease_associations <- disease_associations[n_cases >= 40,]

#### prepare the blood traits and disease used in other sections ####
blood_traits_prs_plot_paper <- readRDS('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/prs/blood_traits_prs_paper.rds')
blood_traits_prs_all_paper <- readRDS('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/prs/blood_traits_prs_all_paper.rds')

blood_traits_prs_plot_paper <- gsub('-', '_', blood_traits_prs_plot_paper)
blood_traits_prs_plot_paper <- gsub('\\.', '_', blood_traits_prs_plot_paper)
blood_traits_prs_all_paper <- gsub('-', '_', blood_traits_prs_all_paper)
blood_traits_prs_all_paper <- gsub('\\.', '_', blood_traits_prs_all_paper)

## use subset of disease associations also used in PRS portion
disease_associations$clinical_trait <- gsub('_', ' ', disease_associations$clinical_trait)
disease_associations <- disease_associations[clinical_trait %in% clinical_annotation_df$name]

# rename the lab associations to be more readable
lab_associations <- merge(lab_associations, lab_trait_annotations[, c('name', 'name2')], by.x = 'clinical_trait', by.y = 'name2')
lab_associations$clinical_trait <- lab_associations$name

lab_associations$name <- NULL
joint_associations <- rbind(lab_associations, disease_associations)
joint_associations$sysmex_pheno <- gsub('-q', '', joint_associations$sysmex_pheno)
joint_associations <- joint_associations[sysmex_pheno %in% blood_traits_prs_all_paper]

## make the sysmex phenotype and clinical trait more readable
joint_associations$sysmex_pheno <- gsub('_', ' ', joint_associations$sysmex_pheno)

# calculate FDR-adjusted q-values the clinical outcomes traits (both disease and labs)
joint_associations
joint_associations$q_sysmex <- qvalue(joint_associations$p_sysmex)$qvalues

hist(log(abs(joint_associations$beta_sysmex)))
hist(log(abs(joint_associations$t_sysmex)))

fwrite(joint_associations, '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/sysmex_readout_clinical_associations_full_2307.tsv', sep='\t')
fwrite(joint_associations[q_sysmex < 0.1,], '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/sysmex_readout_clinical_associations_2307.tsv', sep='\t')

# extract a few columns for plotting
plot_data <- joint_associations[, c('sysmex_pheno', 'clinical_trait',
                                    'p_sysmex', 'q_sysmex', 't_sysmex',
                                    'beta_sysmex', 'channel', 'ptb_name'
)]

#### plot of the full dataset for Supplementary Fig. 2 ####

# reshape into matrix shape
tscore_df <- dcast(plot_data, channel + ptb_name + sysmex_pheno ~ clinical_trait, value.var = 't_sysmex')
tscore_mat <- data.matrix(tscore_df[, -c(1, 2, 3)])
rownames(tscore_mat) <- tscore_df$sysmex_pheno

beta_df <- dcast(plot_data, channel + ptb_name + sysmex_pheno ~ clinical_trait, value.var = 'beta_sysmex')
beta_mat <- data.matrix(beta_df[, -c(1, 2, 3)])
rownames(beta_mat) <- beta_df$sysmex_pheno

pval_df <- dcast(plot_data, channel + ptb_name + sysmex_pheno ~ clinical_trait, value.var = 'q_sysmex')
pval_mat <- data.matrix(pval_df[,-c(1, 2, 3)])
rownames(pval_mat) <- tscore_df$sysmex_pheno

# prepare ordering of rows and columns
## t-score plot

# only showing complete cases
pval_mat <- pval_mat[complete.cases(tscore_mat),]
beta_mat <- beta_mat[complete.cases(beta_mat),]
tscore_mat <- tscore_mat[complete.cases(tscore_mat),]
clinical_annotation_df <- clinical_annotation_df[, c('category', 'type')]

# blood trait colors
bloodtrait_annotation_df <- plot_data[, c('sysmex_pheno', 'channel', 'ptb_name')]
bloodtrait_annotation_df <- as.data.frame(unique(bloodtrait_annotation_df))
rownames(bloodtrait_annotation_df) <- bloodtrait_annotation_df$sysmex_pheno

bloodtrait_annotation_df$sysmex_pheno <- 
  gsub('Baseline', 'Baseline 0h', bloodtrait_annotation_df$sysmex_pheno)
bloodtrait_annotation_df$sysmex_pheno <- 
  gsub('overnight ', '', bloodtrait_annotation_df$sysmex_pheno)

bloodtrait_annotation_df$celltype <- gsub('.*[0-9]h ', '',  bloodtrait_annotation_df$sysmex_pheno)
bloodtrait_annotation_df$celltype <- gsub('[0-9]', '', str_split_fixed(bloodtrait_annotation_df$celltype, ' ', 2)[,1])
bloodtrait_annotation_df <- bloodtrait_annotation_df[, c('channel', 'celltype')]

# define colors for the annotation
ann_colors <- readRDS('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/ann_colors_prs.rds')

# add a single color for hepatic lab values
expanded_category <- c(ann_colors$category, "#5D576B")
names(expanded_category) <- c(names(ann_colors$category), 'Hepatic')
ann_colors$category <- expanded_category

type_cols <- c('#999999', '#a65628')
names(type_cols) <- c('diagnosis', 'lab')
ann_colors[['type']] <- type_cols

callback = function(hc, ...){dendsort(hc, isReverse = TRUE)}
max(abs(tscore_mat))

# make dot font size larger
rg <- 15
p <- pheatmap(
  tscore_mat,
  color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
  breaks = seq(-rg, rg, length.out = 100),
  display_numbers = matrix(ifelse(pval_mat < 0.01, "·", ""), nrow(pval_mat)),
  annotation_row = bloodtrait_annotation_df,
  annotation_col = clinical_annotation_df,
  clustering_distance_cols='correlation',
  angle_col = "45",
  border_color=NA,
  annotation_colors = ann_colors,
  fontsize_row = 3,
  fontsize_number = 10,
  show_colnames=TRUE,
  cutree_rows = 6,
  cutree_cols = 5,
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/sysmex_clinical_associations_2307_full.pdf',
  width = 12,
  height = 16,
  clustering_callback = callback
)
p

qval_labels <- matrix(nrow=nrow(pval_mat), ncol=ncol(pval_mat))
qval_labels[pval_mat <= 1] <- ''
qval_labels[pval_mat <= 0.001] <- '·'
qval_labels[pval_mat <= 0.0001] <- '··'
qval_labels[pval_mat <= 0.00001] <- '···'

colnames(qval_labels) <- colnames(pval_mat)
rownames(qval_labels) <- rownames(pval_mat)

rg <- 1
p <- pheatmap(
  beta_mat,
  color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
  breaks = seq(-rg, rg, length.out = 100),
  display_numbers =qval_labels,
  annotation_row = bloodtrait_annotation_df,
  annotation_col = clinical_annotation_df,
  clustering_distance_cols='correlation',
  angle_col = "45",
  border_color=NA,
  annotation_colors = ann_colors,
  fontsize_row = 3,
  fontsize_number = 10,
  show_colnames=TRUE,
  cutree_rows = 6,
  cutree_cols = 5,
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/sysmex_clinical_associations_2307_full_beta.pdf',
  width = 12,
  height = 16,
  clustering_callback = callback
)
p

#### plot subset of clinical outcomes and blood traits for Figure 3B ####

clinical_trait_counts <- table(joint_associations[q_sysmex < 0.0001]$clinical_trait)
sysmex_pheno_counts <- table(joint_associations[q_sysmex < 0.0001]$sysmex_pheno)

clinical_trait_counts <- clinical_trait_counts[clinical_trait_counts > 2]
sysmex_pheno_counts <- sysmex_pheno_counts[sysmex_pheno_counts > 4]

names(clinical_trait_counts)
names(sysmex_pheno_counts)

tscore_mat_subset <- tscore_mat[names(sysmex_pheno_counts), names(clinical_trait_counts)]
pval_mat_subset <- pval_mat[names(sysmex_pheno_counts), names(clinical_trait_counts)]
beta_mat_subset <- beta_mat[names(sysmex_pheno_counts), names(clinical_trait_counts)]

rownames(tscore_mat_subset) == rownames(pval_mat_subset)
colnames(tscore_mat_subset) == colnames(pval_mat_subset)

sum(is.na(pval_mat_subset))

qval_labels <- matrix(nrow=nrow(pval_mat_subset), ncol=ncol(pval_mat_subset))
qval_labels[pval_mat_subset <= 1] <- ''
qval_labels[pval_mat_subset <= 0.001] <- '·'
qval_labels[pval_mat_subset <= 0.0001] <- '··'
qval_labels[pval_mat_subset <= 0.00001] <- '···'

colnames(qval_labels) <- colnames(pval_mat_subset)
rownames(qval_labels) <- rownames(pval_mat_subset)

# make dot font size larger
rg <- 15
p <- pheatmap(
  tscore_mat_subset,
  color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
  breaks = seq(-rg, rg, length.out = 100),
  display_numbers = qval_labels,
  annotation_row = bloodtrait_annotation_df,
  annotation_col = clinical_annotation_df,
  clustering_distance_cols='correlation',
  angle_col = "45",
  annotation_colors = ann_colors,
  fontsize_row = 5,
  fontsize_col = 6,
  fontsize_number = 8,
  border_color = NA,
  show_colnames=TRUE,
  cutree_rows = 7,
  cutree_cols = 5,
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/sysmex_clinical_associations_2307_subset.pdf',
  width = 8.5,
  height = 8,
  clustering_callback = callback
)
p

rg <- 1
p <- pheatmap(
  beta_mat_subset,
  color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
  breaks = seq(-rg, rg, length.out = 100),
  display_numbers = qval_labels,
  annotation_row = bloodtrait_annotation_df,
  annotation_col = clinical_annotation_df,
  clustering_distance_cols='correlation',
  angle_col = "45",
  annotation_colors = ann_colors,
  fontsize_row = 5,
  fontsize_col = 6,
  fontsize_number = 8,
  border_color = NA,
  show_colnames=TRUE,
  cutree_rows = 7,
  cutree_cols = 5,
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/sysmex_clinical_associations_2307_subset_beta.pdf',
  width = 8.5,
  height = 8,
  clustering_callback = callback
)
p

#### ICA plot, Figure 3C

# which blood traits to show
sel_loadings <- c(
  "wnr Water 15h WBC2 Med FSC",
  "wnr KCl 17h BASO Med SSC",
  "pltf LPS 18h PLT F Med SFL",
  "pltf Pam3CSK4 19h PLT F Percentage",
  "wdf Alhydrogel 21h NE2 NE4 ratio"
)

tscores <- tscore_mat
tscores[is.na(tscores)] <- 0

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

# create a data frame with the top traits
df_traits <- data.frame(IC1 = mixing_matrix[, 1], IC2 = mixing_matrix[, 2], Trait = rownames(mixing_matrix))
df_traits <- df_traits[df_traits$Trait %in% sel_loadings,]

clinical_annotation_df$name <- rownames(clinical_annotation_df)
df <- merge(df, clinical_annotation_df, by.x='Disease', by.y='name')
dark_category_cols <- darken(ann_colors$category, 0.25)
names(dark_category_cols) <- names(ann_colors$category)

# plot the independent components
p <- ggplot(df, aes(x = IC1, y = IC2, color=category)) +
  geom_point() +
  geom_text_repel(aes(label = Disease), segment.colour = NA) +
  labs(x = "IC1", y = "IC2") +
  geom_text_repel(data = df_traits, aes(x = IC1, y = IC2, label = Trait), color = 'black') +
  geom_segment(data = df_traits, 
               aes(x = 0, y = 0, xend = IC1, yend = IC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 'black', 
               alpha = 0.5) +
  scale_color_manual(values=dark_category_cols) +
  theme_minimal()

p
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/raw_readouts_230712_ica_full.pdf',
       width = 10, height = 10)

### subset the disease and blood traits prior to projection
sel_clin_traits <- c(
  unique(disease_associations$clinical_trait),
  'Glucose', 'Glomerular filtration rate', 'QTC interval', 'Creatinine'
)

tscore_df_ica <- dcast(plot_data[clinical_trait %in% sel_clin_traits], channel + ptb_name + sysmex_pheno ~ clinical_trait, value.var = 't_sysmex')
tscore_mat_ica <- data.matrix(tscore_df_ica[, -c(1, 2, 3)])
rownames(tscore_mat_ica) <- tscore_df_ica$sysmex_pheno

tscores <- tscore_mat_ica
tscores[is.na(tscores)] <- 0

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

# create a data frame with the top traits
df_traits <- data.frame(IC1 = mixing_matrix[, 1], IC2 = mixing_matrix[, 2], Trait = rownames(mixing_matrix))

sel_loadings_prs <- c(
  "wnr Water 15h WBC2 Med FSC",
  "wnr KCl 17h BASO Med SSC",
  "pltf LPS 18h PLT F Med SFL",
  "pltf Pam3CSK4 19h PLT F Percentage",
  "wdf Alhydrogel 21h NE2 NE4 ratio"
)

sel_loadings_raw <- c(
  "pltf LPS 18h PLT F Percentage",
  "ret Pam3CSK4 19h RBC SD FSC",
  "ret Captopril 5 5h RET2 Percentage",
  "wdf Alhydrogel 21h NE2 NE4 ratio",
  "wnr KCl 17h WBC Med SSC",
  "pltf Colchicine 20h IPF Med FSC",
  "ret LPS 18h RBC Percentage"
)


heatmap_sysmex_traits <- names(sysmex_pheno_counts)

# align axes to match PRS figure
df_aligned <- df
df_traits_aligned <- df_traits

df_aligned$IC2 <- -df_aligned$IC2
df_traits_aligned$IC2 <- -df_traits_aligned$IC2

df_aligned <- df_aligned[(sqrt(df_aligned$IC1^2 + df_aligned$IC2^2) >= 0.5),]
df_aligned$name <- gsub('_', ' ', df_aligned$Disease)
df_aligned2 <- merge(df_aligned,
                     clinical_annotation_df,
                     by.x='name',
                     by.y='name',
                     how='left')

# plot the independent components
p <- ggplot(df_aligned2, aes(x = IC1, y = IC2, color=category)) +
  geom_point() +
  geom_text_repel(aes(label = Disease), segment.colour = NA) +
  labs(x = "IC1", y = "IC2") +
  geom_text_repel(data = df_traits_aligned[df_traits_aligned$Trait %in% heatmap_sysmex_traits,], aes(x = IC1, y = IC2, label = Trait), color = 'black') +
  geom_segment(data = df_traits_aligned[df_traits_aligned$Trait %in% heatmap_sysmex_traits,], 
               aes(x = 0, y = 0, xend = IC1, yend = IC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 'black', 
               alpha = 0.5) +
  scale_color_manual(values=dark_category_cols) +
  theme_minimal()

p
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/raw_readouts_diagnoses_230712_full_ica.pdf',
       width = 10, height = 10)

p <- ggplot(df_aligned2, aes(x = IC1, y = IC2, color=category)) +
  geom_point() +
  geom_text_repel(aes(label = Disease), segment.colour = NA) +
  labs(x = "IC1", y = "IC2") +
  geom_text_repel(data = df_traits_aligned[df_traits_aligned$Trait %in% sel_loadings_prs,], aes(x = IC1, y = IC2, label = Trait), color = 'black') +
  geom_segment(data = df_traits_aligned[df_traits_aligned$Trait %in% sel_loadings_prs,], 
               aes(x = 0, y = 0, xend = IC1, yend = IC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 'black', 
               alpha = 0.5) +
  scale_color_manual(values=dark_category_cols) +
  theme_minimal()

p
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/raw_readouts_diagnoses_230712_prstraits_ica.pdf',
       width = 10, height = 10)


p <- ggplot(df_aligned2, aes(x = IC1, y = IC2, color=category)) +
  geom_point() +
  geom_text_repel(aes(label = Disease), segment.colour = NA) +
  labs(x = "IC1", y = "IC2") +
  geom_text_repel(data = df_traits_aligned[df_traits_aligned$Trait %in% sel_loadings_raw,], aes(x = IC1, y = IC2, label = Trait), color = 'black') +
  geom_segment(data = df_traits_aligned[df_traits_aligned$Trait %in% sel_loadings_raw,], 
               aes(x = 0, y = 0, xend = IC1, yend = IC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 'black', 
               alpha = 0.5) +
  scale_color_manual(values=dark_category_cols) +
  theme_minimal()

p
ggsave('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/raw_readouts_diagnoses_230712_rawtraits_ica.pdf',
       width = 10, height = 10)

#### show distributions of subset of raw readout - Figure 3A ####

## read in raw readouts and disease/lab status
bmp <- arrow::read_parquet("/mnt/obi0/phi/ehr/obi_biobank_annotations/labs_bmp_21-07-06.parquet")
cmp <- arrow::read_parquet("/mnt/obi0/phi/ehr/obi_biobank_annotations/labs_cmp_21-07-06.parquet")
ecg <- arrow::read_parquet("/mnt/obi0/phi/ehr/obi_biobank_annotations/labs_ecg_21-07-06.parquet")
lipids <- arrow::read_parquet("/mnt/obi0/phi/ehr/obi_biobank_annotations/labs_lp_21-07-06.parquet")
transplant_patients <- arrow::read_parquet('/mnt/obi0/phi/ehr/obi_biobank_annotations/transplant_with_lvad_21-07-06.parquet')

sysmex_readouts <- arrow::read_parquet('/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/all_pheno_projected_filtered_median_wide.parquet')
disease_status <- fread('/mnt/obi0/phi/ehr/obi_biobank_annotations/disease_230712.csv', sep='\t')
covariates <-  arrow::read_parquet("/mnt/obi0/phi/ehr/obi_biobank_annotations/covariates_21-07-06.parquet")

outcomes_df <- Reduce(function(x, y) 
  merge(x, y, by = "PatientID", all = TRUE),
  list(bmp, cmp, ecg, lipids, transplant_patients, disease_status, covariates)
)

plot_readout_outcomes <- function(
    readout,
    outcome,
    name,
    readout_label,
    ymin = -Inf,
    ymax = Inf
) {
  
  perturb_readouts <- sysmex_readouts[, c('PatientID', readout)]
  perturb_readouts <- na.exclude(perturb_readouts)
  perturb_outcomes <- merge(perturb_readouts, outcomes_df)
  perturb_outcomes$outcome <- perturb_outcomes[[outcome]]
  perturb_outcomes$readout <- perturb_outcomes[[readout]]
  perturb_outcomes$age <- as.numeric(perturb_outcomes$rc_consent_age)
  
  covariates <- c('SexDSC', 'age', 'Race')
  
  plot_outcomes <- na.exclude(perturb_outcomes[, c('readout', 'outcome', covariates)])
  
  # filter outliers / measurement errors if ymin/ymax are supplied
  plot_outcomes <- plot_outcomes[plot_outcomes$readout > ymin,]
  plot_outcomes <- plot_outcomes[plot_outcomes$readout < ymax,]
  
  p <- ggplot(plot_outcomes, aes(x=readout, y=outcome, color=age, fill=age, group=outcome)) +
    stat_slab(alpha = .5, height=1.05) +
    geom_dots(smooth = "unbounded", position = "dodgejust") +
    stat_pointinterval(position = position_dodge(width = .1, preserve = "single")) +
    facet_grid('SexDSC~.') +
    theme_minimal() +
    ggtitle(paste0(name, ' (n=', sum(plot_outcomes$outcome == TRUE), '/', length(plot_outcomes$outcome), ')')) +
    xlab(readout_label) +
    coord_flip()
  p
}

p <-
plot_readout_outcomes(
  readout = "('LPS 18h', 'ret_RBC_Percentage')",
  outcome = 'purpura and haemorrhagic conditions',
  name = 'Purpura and haemorrhagic conditions',
  readout_label = 'ret LPS RBC Percentage',
  ymin = 85
  )
ggsave(
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/readout_associations/ret_lps_rbc_purpura.pdf', 
  plot=p, width = 6, height = 5
)

readout_names <- 
colnames(sysmex_readouts)

p <-
  plot_readout_outcomes(
    readout = "('Alhydrogel 21h', 'wdf_NE2_NE4_ratio')",
    outcome = 'type 2 diabetes mellitus',
    name = 'Type 2 diabetes mellitus',
    readout_label = 'wdf Alhydrogel NE2-NE4 ratio',
    #ymin = 85
    #ymax=0.6
  )
p
ggsave(
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/readout_associations/wdf_alhydrogel_ne2_ne4_t2d.pdf', 
  plot=p, width = 6, height = 5
)

p <-
  plot_readout_outcomes(
    readout = "('Alhydrogel 21h', 'wdf_NE4_SD_SFL')",
    outcome = 'type 2 diabetes mellitus',
    name = 'Type 2 diabetes mellitus',
    readout_label = 'wdf Alhydrogel NE4 SD SFL',
    #ymin = 85
    #ymax=0.6
  )
p
ggsave(
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/readout_associations/wdf_alhydrogel_ne4_sd_sfl_t2d.pdf', 
  plot=p, width = 6, height = 5
)

p <-
  plot_readout_outcomes(
    readout = "('KCl 17h', 'wnr_WBC_Med_SSC')",
    outcome = 'chronic kidney disease',
    name = 'chronic kidney disease',
    readout_label = 'wnr KCl WBC Med SSC',
    #ymin = 85
    #ymax=0.6
  )
p
ggsave(
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/readout_associations/wnr_kcl_WBC_SSC_ckd.pdf', 
  plot=p, width = 6, height = 5
)

p <-
  plot_readout_outcomes(
    readout = "('Pam3CSK4 19h', 'ret_RBC_SD_FSC')",
    outcome = 'heart failure',
    name = 'Heart failure',
    readout_label = 'ret Pam3CSK4 RBC SD FSC',
    #ymin = 85
    #ymax=28
  )
p
ggsave(
  filename = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/plots/readout_associations/ret_RBC_SD_FSC_heart_failure.pdf', 
  plot=p, width = 6, height = 5
)
