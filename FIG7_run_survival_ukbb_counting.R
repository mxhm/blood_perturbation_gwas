# Figure 7A, 7B: run survival analyses using Cox PH models with delayed entry (counting process-style input)
# Use blood-PGS as part of the models, and clinical outcomes based on date of first diagnosis

library(data.table)
library(ukbtools)
library(matrixStats)
library(lubridate)
library(parallel)
library(stringr)

library(ggplot2)
library(survminer)
library(survival)

# Rscript /home/mhomilius/projects/sysmex/scripts/r_plotting_revision/fig7_run_survival_ukbb_counting_230706.R

covars <- c(
  'eid',
  "sex_f31_0_0", 
  "year_of_birth_f34_0_0",
  "ethnic_background_f21000_0_0",
  "ethnic_background_f21000_1_0",
  "ethnic_background_f21000_2_0",
  "age_2018", "age_2023",
  "birth_year_date",
  "age_at_recruitment_f21022_0_0",
  'date_of_attending_assessment_centre_f53_0_0',
  'date_lost_to_followup_f191_0_0',
  'age_entry_assessment',
  'age_lost_followup',
  'age_first_icd',
  'earliest_icd_date',
  'age_entry'
)


# self-reported age for a small number of traits
phenotypes_age <- c(
  'age_diabetes_diagnosed_f2976',
  'age_asthma_diagnosed_f3786',
  'age_hay_fever_rhinitis_or_eczema_diagnosed_f3761',
  'age_angina_diagnosed_f3627',
  'age_deepvein_thrombosis_dvt_blood_clot_in_leg_diagnosed_f4012',
  'age_emphysemachronic_bronchitis_diagnosed_f3992',
  'age_heart_attack_diagnosed_f3894',
  'age_pulmonary_embolism_blood_clot_in_lung_diagnosed_f4022',
  'age_stroke_diagnosed_f4056'
)

# note: list below contains clinical outcomes not used in paper
# first occurence in multiple sources
phenotypes_icd_date <- c(
  'date_d50_first_reported_iron_deficiency_anaemia_f130622_0_0',
  'date_d59_first_reported_acquired_haemolytic_anaemia_f130638_0_0',
  'date_d63_first_reported_anaemia_in_chronic_diseases_classified_elsewhere_f130646_0_0',
  'date_d69_first_reported_purpura_and_other_haemorrhagic_conditions_f130658_0_0',
  
  'date_e10_first_reported_insulindependent_diabetes_mellitus_f130706_0_0',
  'date_e11_first_reported_noninsulindependent_diabetes_mellitus_f130708_0_0',
  'date_e66_first_reported_obesity_f130792_0_0',
  'date_e75_first_reported_disorders_of_sphingolipid_metabolism_and_other_lipid_storage_disorders_f130808_0_0',
  "date_e78_first_reported_disorders_of_lipoprotein_metabolism_and_other_lipidaemias_f130814_0_0",
  
  "date_f01_first_reported_vascular_dementia_f130838_0_0",
  "date_f03_first_reported_unspecified_dementia_f130842_0_0",
  "date_f40_first_reported_phobic_anxiety_disorders_f130904_0_0",
  "date_f41_first_reported_other_anxiety_disorders_f130906_0_0",

  'date_g31_first_reported_other_degenerative_diseases_of_nervous_system_not_elsewhere_classified_f131038_0_0',
  "date_g43_first_reported_migraine_f131052_0_0",
  'date_g45_first_reported_transient_cerebral_ischaemic_attacks_and_related_syndromes_f131056_0_0',
  'date_g56_first_reported_mononeuropathies_of_upper_limb_f131074_0_0',
  'date_g58_first_reported_other_mononeuropathies_f131078_0_0',

  'date_h81_first_reported_disorders_of_vestibular_function_f131252_0_0',
  'date_h34_first_reported_retinal_vascular_occlusions_f131180_0_0',
  
  "date_i10_first_reported_essential_primary_hypertension_f131286_0_0",
  "date_i20_first_reported_angina_pectoris_f131296_0_0",
  "date_i21_first_reported_acute_myocardial_infarction_f131298_0_0",
  "date_i25_first_reported_chronic_ischaemic_heart_disease_f131306_0_0",
  'date_i35_first_reported_nonrheumatic_aortic_valve_disorders_f131324_0_0',
  'date_i40_first_reported_acute_myocarditis_f131334_0_0',
  'date_i42_first_reported_cardiomyopathy_f131338_0_0',
  'date_i48_first_reported_atrial_fibrillation_and_flutter_f131350_0_0',
  "date_i50_first_reported_heart_failure_f131354_0_0",
  'date_i63_first_reported_cerebral_infarction_f131366_0_0',
  "date_i70_first_reported_atherosclerosis_f131380_0_0",
  'date_i71_first_reported_aortic_aneurysm_and_dissection_f131382_0_0',
  "date_i77_first_reported_other_disorders_of_arteries_and_arterioles_f131390_0_0",
  'date_i82_first_reported_other_venous_embolism_and_thrombosis_f131400_0_0',
  'date_i83_first_reported_varicose_veins_of_lower_extremities_f131402_0_0',
  
  "date_j44_first_reported_other_chronic_obstructive_pulmonary_disease_f131492_0_0",
  "date_j45_first_reported_asthma_f131494_0_0",
  "date_j46_first_reported_status_asthmaticus_f131496_0_0",
  
  'date_k52_first_reported_other_noninfective_gastroenteritis_and_colitis_f131630_0_0',
  "date_k58_first_reported_irritable_bowel_syndrome_f131638_0_0",
  'date_k85_first_reported_acute_pancreatitis_f131682_0_0',
  
  'date_m08_first_reported_juvenile_arthritis_f131854_0_0',
  'date_m10_first_reported_gout_f131858_0_0',
  "date_m31_first_reported_other_necrotising_vasculopathies_f131892_0_0",
  "date_m32_first_reported_systemic_lupus_erythematosus_f131894_0_0",
  'date_m34_first_reported_systemic_sclerosis_f131898_0_0',
  "date_m72_first_reported_fibroblastic_disorders_f131950_0_0",
  'date_m79_first_reported_other_soft_tissue_disorders_not_elsewhere_classified_f131960_0_0',
  
  'date_n10_first_reported_acute_tubulointerstitial_nephritis_f132016_0_0',
  'date_n12_first_reported_tubulointerstitial_nephritis_not_specified_as_acute_or_chronic_f132020_0_0',
  'date_n14_first_reported_drug_and_heavymetalinduced_tubulointerstitial_and_tubular_conditions_f132024_0_0',
  "date_n17_first_reported_acute_renal_failure_f132030_0_0",
  "date_n18_first_reported_chronic_renal_failure_f132032_0_0",
  "date_n19_first_reported_unspecified_renal_failure_f132034_0_0",
 
  'date_l52_first_reported_erythema_nodosum_f131758_0_0',
  "date_l93_first_reported_lupus_erythematosus_f131828_0_0",
  
  'date_o24_first_reported_diabetes_mellitus_in_pregnancy_f132202_0_0'
)

phenotypes_algo_date <- c(
  'date_of_asthma_report_f42014_0_0',
  'date_of_all_cause_dementia_report_f42018_0_0',
  'date_of_vascular_dementia_report_f42022_0_0',
  "date_of_frontotemporal_dementia_report_f42024_0_0",
  'date_of_chronic_obstructive_pulmonary_disease_report_f42016_0_0',
  'date_of_end_stage_renal_disease_report_f42026_0_0',
  'date_of_myocardial_infarction_f42000_0_0',
  'date_of_stemi_f42002_0_0',
  'date_of_nstemi_f42004_0_0',
  'date_of_stroke_f42006_0_0',
  'date_of_ischaemic_stroke_f42008_0_0',
  'date_of_intracerebral_haemorrhage_f42010_0_0',
  'date_of_subarachnoid_haemorrhage_f42012_0_0'
)

phenotypes_age_computed <- pheno_age <- gsub('date_', 'age_', c(phenotypes_icd_date, phenotypes_algo_date))
phenotypes <- c(phenotypes_age, phenotypes_age_computed)
covar_pheno <- c(covars, phenotypes_icd_date, phenotypes_algo_date, phenotypes)

# #### preprocessing - run only once ####
# ukb_data_metab <- ukb_df("ukb672103", path = '/mnt/obi0/mhomilius/projects/sysmex/data/UKBB')
# ukb_data_diag <- ukb_df("ukb671218", path = '/mnt/obi0/mhomilius/projects/sysmex/data/UKBB')
# ukb_data_diag_dt <- data.table(ukb_data_diag)
# 
# # # icd coding file is here:
# # ukbb_icd10_mapping <- fread('/mnt/obi0/mhomilius/projects/sysmex/data/UKBB/ukbb_icd10_mapping.tsv')
# # sysmex_icd10_codes <- readxl::read_excel('/mnt/obi0/mhomilius/projects/sysmex/data/atlas/clinical_traits_icd10.xlsx')
# 
# # add agreggate pheno values across the different visits for self-reported traits
# for (pheno in phenotypes_age){
#   print(pheno)
#   pheno_visits <- paste0(pheno, c('_0_0', '_1_0', '_2_0', '_3_0'))
#   ukb_data_diag_dt[[pheno]] <- ukb_data_diag_dt[,..pheno_visits][, (pheno) := rowMedians(as.matrix(.SD), na.rm = TRUE)][[pheno]]
# }
# 
# # calculcate earliest icd10 date across all ICD10 codes
# icd10_cols <- all_cols[grep(pattern = 'date_.[0-9]*_first_reported', x = all_cols)]
# select_cols <- c('eid', icd10_cols)
# icd_dates <- ukb_data_diag_dt[, ..select_cols]
# 
# # convert everything to dates
# icd_dates[ , (icd10_cols) := lapply(.SD, as.Date), .SDcols = icd10_cols]
# 
# # find the earliest ICD10 report date per subject
# icd_dates$earliest_icd_date <- apply(icd_dates[, .SD, .SDcols = icd10_cols], 1, min, na.rm = TRUE)
# 
# ukb_data_diag_dt$earliest_icd_date <- icd_dates$earliest_icd_date
# # calculate age at diagnosis for date-based traits
# 
# # use age in 2023 as end date for controls
# ukb_data_diag_dt$birth_year_date <- as_date(as.character(ukb_data_diag_dt$year_of_birth_f34_0_0), format='%Y')
# ukb_data_diag_dt$age_2018 <- as.numeric(interval(ukb_data_diag_dt$birth_year_date, as_date('2018-01-01')), 'years')
# ukb_data_diag_dt$age_2023 <- as.numeric(interval(ukb_data_diag_dt$birth_year_date, as_date('2023-01-01')), 'years')
# ukb_data_diag_dt$age_entry_assessment <- as.numeric(interval(
#   ukb_data_diag_dt$birth_year_date,
#   as_date(ukb_data_diag_dt$date_of_attending_assessment_centre_f53_0_0)), 'years')
# 
# ukb_data_diag_dt$age_first_icd <- as.numeric(interval(
#   ukb_data_diag_dt$birth_year_date,
#   ukb_data_diag_dt$earliest_icd_date), 'years')
# 
# ukb_data_diag_dt$age_lost_followup <- as.numeric(interval(
#   ukb_data_diag_dt$birth_year_date,
#   as_date(ukb_data_diag_dt$date_lost_to_followup_f191_0_0)), 'years')
# 
# # use the earliest ICD10 date or the assessment center date as entry date, whichever is earlier
# ukb_data_diag_dt$age_entry <- pmin(ukb_data_diag_dt$age_first_icd, ukb_data_diag_dt$age_entry_assessment, na.rm = TRUE)
# 
# #
# #hist(ukb_data_diag_dt$age_first_icd)
# #hist(ukb_data_diag_dt$age_entry_assessment)
# #hist(ukb_data_diag_dt$age_entry)
# 
# # add age-based values for date-based traits
# for (pheno in c(phenotypes_icd_date, phenotypes_algo_date)){
#   # TODO need to remove subjects with special values when running survival models
#   pheno_age <- gsub('date_', 'age_', pheno)
#   ukb_data_diag_dt[[pheno_age]] <- as.numeric(interval(ukb_data_diag_dt$birth_year_date, ukb_data_diag_dt[[pheno]]), 'years')
# }
# 
# # create a dataframe with subset of data needed for the survival models
# pheno_subset <- ukb_data_diag_dt[, ..covar_pheno]
# pheno_subset$eid <- as.character(pheno_subset$eid)
# 
# fwrite(pheno_subset, file='/mnt/obi0/mhomilius/projects/sysmex/data/UKBB/pheno_covariate_subset_expanded_2307.tsv', sep='\t')

#### survival analyses ####
pheno_subset <- fread('/mnt/obi0/mhomilius/projects/sysmex/data/UKBB/pheno_covariate_subset_expanded_2307.tsv')
pheno_subset$eid <- as.character(pheno_subset$eid)

# # writing out subject ids of white/eur ancestry subjects for PCA filtering
# extract_cols <- c('eid', 'eid')
# fwrite(
#   pheno_subset[ethnic_background_f21000_0_0 == '1001', ..extract_cols],
#   file='~/projects/sysmex/data/UKBB/genotyped/eur_subjects',
#   col.names = FALSE,
#   row.names = FALSE,
#   sep='\t'
#   )

pca_covariates <- fread('/mnt/obi0/mhomilius/projects/sysmex/data/UKBB/genotyped/ukb_cal_allChrs_qc.pruned.eur.pca.eigenvec')
pca_covariates$eid <- as.character(pca_covariates$IID)
pheno_subset <- merge(pheno_subset, pca_covariates, by.x = 'eid', by.y='eid')

# pca_covariates[1:10000,]
# library(ggplot2)
# ggplot(pca_covariates[1:10000,], aes(x=PC1, y=PC2)) +
#   geom_point()
# 
# ggplot(pca_covariates[1:10000,], aes(x=PC3, y=PC4)) +
#   geom_point()

plotKM <-
  function(fit,
           data,
           disease_label,
           output_file,
           start_year = 0,
           end_year = 15) {
    p <- ggsurvplot(
      fit,
      data = data,
      risk.table = TRUE,
      pval = TRUE,
      pval.coord = c(start_year, 0.05),
      risk.table.height = 0.25,
      xlim = c(start_year - 0.1, end_year + 0.1),
      break.time.by = 10,
      #palette = c('plum3', 'skyblue3', 'darkseagreen', 'firebrick', 'goldenrod3'),
      risk.table.pos = 'out',
      test.for.trend = TRUE,
      censor = FALSE
    )
    p$plot <- p$plot +
      xlab('Time (years)') +
      theme_bw()
    print(p$plot)
    ggsave(paste0(output_file, '_tte.pdf'),
           width = 8,
           height = 8)
    
    p$table <- p$table + xlab("Time (years)")
    print(p$table)
    ggsave(paste0(output_file, '_table.pdf'),
           width = 8,
           height = 2)
    
    p <- ggsurvplot(
      fit,
      data = data,
      risk.table = FALSE,
      pval = TRUE,
      pval.coord = c(start_year, 0.755),
      xlim = c(start_year - 0.1, end_year + 0.1),
      ylim = c(0.75, 1),
      break.time.by = 10,
      #palette = c('plum3', 'skyblue3', 'darkseagreen', 'firebrick', 'goldenrod3'),
      risk.table.pos = 'out',
      test.for.trend = TRUE,
      # CHANGED testing for trend among groups
      censor = FALSE
    )
    p$plot <- p$plot +
      xlab('Time (years)') +
      theme_bw()
    plot(p$plot)
    ggsave(paste0(output_file, '_inset.pdf'),
           width = 8,
           height = 4)
  }

logrank_tft_zscore <- function(fit, data) {
  # code below extracted from surv_pvalue and survminer:::.pvalue
  method <- "1"
  tenfit <- survMisc::ten(eval(fit$call$formula), data = data)
  null_dev <- capture.output(survMisc::comp(tenfit))
  tests <- attributes(tenfit)$tft$tft
  
  lr_p <- tests$pChisq
  lr_z <- tests$Z
  lr_w <- tests$W
  
  # TODO hardcoded
  list(
    method = 'Log-rank, tft',
    term = 'quartile',
    p.value = lr_p[lr_w == method],
    statistic = lr_z[lr_w == method]
  )
}


overwrite <- FALSE
run_survival_analysis <- function(score_file) {
  output_root <-
    '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/survival_ct_p01_sex_pc2_counting_ukbb230703/'
  blood_trait <-
    str_split_fixed(basename(score_file), '-q', n = 2)[, 1]
  
  if (!file.exists(paste0(output_root, blood_trait, '.stats.txt')) ||
      overwrite) {
    file.create(paste0(output_root, blood_trait, '.stats.txt'))
    
    print(score_file)
    
    
    stats_ <-
      mclapply(phenotypes, mc.cores = 8, function(pheno) {
        
        clinical_trait <- pheno
        output_file <- paste0(blood_trait, '|', clinical_trait)
        output_path <- paste0(output_root, output_file)
        stats_file <- paste0(output_path, '--stats.tsv')
        
        if (!file.exists(stats_file) || overwrite) {
          stats_df <- NULL
          
          # create empty file
          file.create(stats_file)
          print(output_file)
          
          score_df <- fread(score_file)
          score_df$IID <- as.character(score_df$IID)
          
          # scale all the risk scores so that beta values are comparable
          score_df$risk_score <- scale(score_df$SCORE)
          pheno_df <-
            merge(pheno_subset,
                  score_df,
                  by.x = 'eid',
                  by.y = 'IID')
          
          
          # subset to european ancestry cohort
          pheno_df <- pheno_df[ethnic_background_f21000_0_0 == '1001', ]
          pheno_df$Sex <- as.factor(pheno_df$sex_f31_0_0)
          
          pheno_df$case <- !is.na(pheno_df[[pheno]])
          
          # only proceed if at least 40 cases
          if (sum(pheno_df$case) >= 20) {
            disease_label <- ''
            
            pheno_df$years_to_event <- pheno_df$age_2023
            # if lost to follow up, use that end data
            pheno_df$age_lost_followup2 <- as.numeric(interval(
              as_date(pheno_df$birth_year_date), as_date(pheno_df$date_lost_to_followup_f191_0_0)), 'years')
            
            lost_follow_up <- !is.na(pheno_df$age_lost_followup2)
            pheno_df[lost_follow_up,]$years_to_event <- pheno_df[lost_follow_up,]$age_lost_followup2
            
            # for cases, use the diagnosis date
            pheno_df[pheno_df$case, 'years_to_event'] <- pheno_df[pheno_df$case,][[pheno]]
            
            # drop anyone with a special date if this was a date-based trait
            
            # 1900-01-01 date unknown
            
            # ICD10 traits
            # 1901-01-01 before birth
            # 1902-02-02 matching birth
            # 1903-03-03 birth year
            # 2037-07-07 event date in the future
            
            # 1903-03-03
            
            date_pheno <- gsub('age', 'date', pheno)
            if (date_pheno %in% colnames(pheno_df)){
              special_dates <- pheno_df[[date_pheno]] %in% as.Date(c('1900-01-01', '1901-01-01', '1902-02-02', '1903-03-03', '2037-07-07'))
              #special_dates <- pheno_df[[date_pheno]] %in% c('1900-01-01', '1901-01-01', '1902-02-02', '1903-03-03', '2037-07-07')
              pheno_df <- pheno_df[!special_dates]
            }
            
            hist(pheno_df$years_to_event)
            summary(pheno_df$years_to_event)
            

            pheno_df <- pheno_df[years_to_event >= 0,]
            print('Median follow up (years after birth)')
            print(median(pheno_df$years_to_event))
            print(summary(pheno_df$years_to_event))
            
            # create quartiles and fit survival model
            qts = quantile(pheno_df$risk_score)
            print(blood_trait)
            print(qts)
            pheno_df$quartile <-
              with(pheno_df,
                   cut(
                     risk_score,
                     qts,
                     include = TRUE,
                     labels = FALSE
                   ))
            
            #### prepare counting process-style input for coxph
            
            # there are some entries were a special code is used to code the birth year
            pheno_df$age_entry2 <- pmax(pheno_df$age_entry, 0)
            
            # start observation interval
            pheno_df$time <- pheno_df$age_entry2
            pheno_df$status <- 0
            
            # for cases, use time to indicate event time
            pheno_df[case == TRUE]$time <- pheno_df[case == TRUE]$years_to_event
            pheno_df[case == TRUE]$status <- 1
            
            # for cases, where the entry time is less than half a year before event time, set to delayed entry
            pheno_df[((case == TRUE) & (years_to_event - age_entry2 <= 1))]$status <- 2
            
            # all of the observations stop at years_to_event
            pheno_df$time2 <- pheno_df$years_to_event
            
            data <- pheno_df[, c('eid', 'age_entry2', 'years_to_event', 'case', 'Sex', 'PC1', 'PC2', 'PC3', 'risk_score', 'quartile')]
            setnames(data,  c('eid', 'age_entry2', 'years_to_event', 'case'), c('id', 'time1', 'time2', 'status'))
            
            # drop ~1k controls with entry time at birth
            data <- data[!((status == FALSE) & (time1 == 0)),]
            
            # Create a data frame for the period from the start of the study until each individual enters
            data_start <- data.frame(id = data$id, 
                                     start = rep(0, nrow(data)),  # start of study
                                     stop = data$time1,  # when the individual enters
                                     event = rep(0, nrow(data)),  # no event
                                     Sex = data$Sex,
                                     PC1 = data$PC1,
                                     PC2 = data$PC2,
                                     PC3 = data$PC3,
                                     risk_score = data$risk_score, 
                                     quartile=data$quartile)
            
            # Create a data frame for the period from when each individual enters until they have the event or are censored
            data_stop <- data.frame(id = data$id, 
                                    start = data$time1,  # when the individual enters
                                    stop = data$time2,  # when the individual has the event or is censored
                                    event = as.numeric(data$status == 1),  # whether the event happened
                                    Sex = data$Sex,
                                    PC1 = data$PC1,
                                    PC2 = data$PC2,
                                    PC3 = data$PC3,
                                    risk_score = data$risk_score,
                                    quartile=data$quartile)
            
            # Adjust the 'event' values in 'data_start' to handle instant events
            # Include events within 1 year of entry date as instant events
            instant_events <- which((data_stop$start - data_stop$stop) > -1 & data_stop$event == 1)
            data_start$event[instant_events] <- 1  # set event = 1
            
            # for some self-reported traits, the event date is much earlier than the earliest ICD diagnosis
            data_stop$start[instant_events] <- data_stop$stop[instant_events] - 0.01 
            
            # Combine the two data frames
            data_counting <- rbind(data_start, data_stop)
            data_stop$event[instant_events] <- 0  # set event = 0
            data_counting2 <- rbind(data_start, data_stop)
            
            # drop controls where stop and start are 0
            data_counting <- data_counting[!(data_counting$start >= data_counting$stop),]
            data_counting2 <- data_counting2[!(data_counting2$start >= data_counting2$stop),]
            
            # fit the Cox models (using only fit3 for paper, fit1 and fit2 double-count the events)
            #fit1 <- coxph(Surv(start, stop, event, type = "counting") ~ Sex + PC1 + PC2 + risk_score, data = data_counting)
            #fit2 <- coxph(Surv(start, stop, event, type = "counting") ~ strata(Sex) + PC1 + PC2 + risk_score, data = data_counting)
            fit3 <- coxph(Surv(start, stop, event, type = "counting") ~ Sex + PC1 + PC2 + risk_score, data = data_counting2)
            #fit4 <- coxph(Surv(start, stop, event, type = "counting") ~ strata(Sex) + PC1 + PC2 + risk_score, data = data_counting2)
            
            # Print the summary
            # print(summary(fit1))
            # print(summary(fit2))
            # print(summary(fit3))
            # print(summary(fit4))
            
            fit_to_stats <- function(fit, name){
              summ_obj <- summary(fit)
              stats_df <- data.frame(broom::tidy(fit))
              stats_df$method = name
              stats_df$concordance <- summ_obj$concordance['C']
              stats_df$se_concordance <- summ_obj$concordance['se(C)']
              stats_df$case_count <- sum(pheno_df$case)
              stats_df$blood_trait <- blood_trait
              stats_df$clinical_trait <- clinical_trait
              stats_df
            }
            
            #stats_df1 <- fit_to_stats(fit1, 'fit1')
            #stats_df2 <- fit_to_stats(fit2, 'fit2')
            stats_df3 <- fit_to_stats(fit3, 'fit3')
            #stats_df4 <- fit_to_stats(fit4, 'fit4')
            
            #stats_df <- rbindlist(list(stats_df1, stats_df2, stats_df3, stats_df4))
            stats_df <- stats_df3
            pval <- stats_df[(stats_df$term == 'risk_score') & (method == 'fit3'), ]$p.value
            print(stats_df)

            # plot if cox model with covariates is significant
            if (!is.na(pval) & pval < 0.01) {
              
              ggplot(pheno_df, aes(x=years_to_event - age_entry2, fill=as.factor(status), group=status)) +
                geom_histogram(position = 'stack') +
                scale_fill_brewer(palette='Set2')
              ggsave(file=gsub('--stats.tsv', '--entry_event_time.pdf', stats_file), width = 6, height = 4)
              ggplot(pheno_df, aes(x=years_to_event, fill=as.factor(status), group=status)) + geom_histogram() + scale_fill_brewer(palette='Set2')
              ggsave(file=gsub('--stats.tsv', '--status_age.pdf', stats_file), width = 6, height = 4)
              
              # subset cases/controls for plotting if they are very large
              pheno_df_plot <- rbind(pheno_df[(case)][sample(.N, min(.N, 30000))],
                                     pheno_df[!(case)][sample(.N,50000)])
              
              sfit <- survfit(Surv(years_to_event, case) ~ quartile, data = pheno_df_plot)
              
              plotKM(
                sfit,
                data = pheno_df_plot,
                output_file = output_path,
                disease_label = disease_label,
                start_year = 0,
                end_year = 85
              )
            }
            
            fwrite(file = stats_file, stats_df, sep = '\t')
            stats_df
          }
        }
      }
    )
    stats_df <- rbindlist(stats_)
    fwrite(stats_df,
           paste0(output_root, blood_trait, '.stats.txt'),
           sep = '\t')
  }
}

prs_files <- list.files(path = '/mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/prs/plink_prs_ukbb221031/',
                        pattern = '*.joint.prs.0.1.profile',
                        full.names = TRUE)

mclapply(prs_files, mc.cores = 4, run_survival_analysis)


# cd /mnt/obi0/phi/gwas/gwas_analyses/sysmex_custom_gates_v9-obi2020_10_28/revision/
# Rscript /home/mhomilius/projects/sysmex/scripts/r_plotting_revision/prepare_ukbb_phenotypes_covariates_delayed_entry.R
