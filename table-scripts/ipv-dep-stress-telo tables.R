rm(list=ls())

library('flextable')
library('officer')
source(here::here("0-config.R"))
library('here')
library('data.table')

# load enrollment characteristics and results
d <- read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-stress-growth-covariates-stresslab-anthro.csv"))
H1 <- readRDS(here('results/unadjusted/H1_res.RDS'))
H2 <- readRDS(here('results/unadjusted/H2_res.RDS'))
H3 <- readRDS(here('results/unadjusted/H3_res.RDS'))
H1adj <- readRDS(here('results/adjusted/H1_adj_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_res.RDS'))


#### Functions for growth tables ####
source(here::here("table-functions.R"))
  # format for export
  flextbl<-flextable(tbl, col_keys=names(tbl))
  flextbl <- set_header_labels(flextbl,
                               values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                             "V6" = "Predicted Outcome at 25th Percentile", "V7" = "Predicted Outcome at 75th Percentile", "V8" = "Coefficient (95% CI)", "V9" = "P-value", "V10" = "FDR Corrected P-value",
                                             "V11" = "Predicted Outcome at 25th Percentile", "V12" = "Predicted Outcome at 75th Percentile", "V13" = "Coefficient (95% CI)", "V14" = "P-value", "V15" = "FDR Corrected P-value"))
  flextbl <- add_header_row(flextbl, values = c("","","","","", "Unadjusted", "Fully adjusted"), colwidths=c(1,1,1,1,1,5,5))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,10))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
  flextbl <- autofit(flextbl, part = "all")
  flextbl <- fit_to_width(flextbl, max_width=8)
  flextbl
}


#### MAIN TABLES ####
#### Table 1 ####
# Characteristics of participants
nperc <- function(vector){
  n <- sum(vector==1, na.rm=T)
  perc <- round(n/sum(!is.na(vector))*100)
  paste(n, " (", perc, "%)", sep="")
}

mediqr <- function(vector){
  quantiles <- round(quantile(vector, na.rm=T), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

n_med_col <- c(nperc(d$sex), mediqr(d$t2_f2_8ip), mediqr(d$t2_f2_23d), mediqr(d$t2_f2_VI), mediqr(d$t2_f2_12i),
               mediqr(d$t3_cort_slope), mediqr(d$t3_residual_cort), mediqr(d$t3_saa_slope), mediqr(d$t3_residual_saa),
               mediqr(d$t3_map), mediqr(d$t3_hr_mean), mediqr(d$t3_gcr_mean), mediqr(d$t3_gcr_cpg12),
               mediqr(d$laz_t2), mediqr(d$waz_t2), mediqr(d$whz_t2), mediqr(d$hcz_t2),
               mediqr(d$laz_t3), mediqr(d$waz_t3), mediqr(d$whz_t3), mediqr(d$hcz_t3),
               nperc(d$diar7d_t2), nperc(d$diar7d_t3), mediqr(d$momage), mediqr(d$momheight), 
               mediqr(d$momeduy), mediqr(d$cesd_sum_t2), mediqr(d$cesd_sum_ee_t3), mediqr(d$pss_sum_mom_t3), 
               nperc(d$life_viol_any_t3))

tbl1 <- data.table("C1" = c("Child","","","","","","","","","","","","","","","","","","","","","","","Mother","","","","","",""),
                   "C2" = c("", "Urinary F2-isoprostanes (Year 1)","","","", "Salivary cortisol reactivity (Year 2)","", "sAA reactivity (Year 2)","",
                           "SAM biomarkers (Year 2)","", "Glucocorticoid receptor","", "Anthropometry (14 months, Year 1)","","","",
                           "Anthropometry (28 months, Year 2)","","","", "Diarrhea (14 months, Year 1)", "Diarrhea (28 months, Year 2)","",
                           "Anthropometry at enrollment", "Education", "Depression at Year 1", "Depression at Year 2", "Perceived stress at Year 2", 
                           "Intimate partner violence"),
                   "C3" = c("Female", "iPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a-VI", "8,12-iso-iPF(2a)-VI", 
                           "Change in slope between pre- and post-stressor cortisol", "Cortisol residualized gain score", 
                           "Change in slope between pre- and post-stressor sAA change", "sAA residualized gain score",
                           "Mean arterial pressure", "Resting heart rate", "NR3C1 exon 1F promoter methylation", "NGFI-A transcription factor binding site methylation",
                           "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                           "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                           "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
                           "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
                   "C4" = n_med_col)

tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
tbl1flex <- set_header_labels(tbl1flex,
                        values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or median (IQR)"))
tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
tbl1flex <- autofit(tbl1flex, part = "all")
tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
tbl1flex <- fit_to_width(tbl1flex, max_width=8)
names(tbl1)<- c("","","","n (%) or median (IQR)")


#### Table 2 ####

exposure <- c("life_viol_any_t3","viol_any_preg", "viol_any_t2","viol_12m_any_t3")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Lifetime Exposure to IPV Year 2","Exposure to IPV during Pregnancy",
              "Exposure to IPV during first year of child's life", "Exposure to IPV in Past 12 Months Year 2")
out_var <- c("Telomere length Z-score Year 1", "Telomere length Z-score Year 2","Change in Telomere length Z-Score")

tbl2 <- growth_tbl("IPV", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl2flex <- growth_tbl_flex("IPV", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl2supp <- growth_tbl("IPV", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl2flexsupp <- growth_tbl_flex("IPV", expo_var, out_var, exposure, outcome, H1, H1adj)

#### Table 3 ####
exposure <- c("cesd_sum_t2", "cesd_sum_t2_binary", "cesd_sum_ee_t3", "cesd_sum_ee_t3_binary")
outcome <- c("TS_t2_Z","TS_t3_Z","delta_TS_Z")
expo_var <- c("Maternal Depression Year 1","Binary Maternal Depression Year 1",
              "Maternal Depression Year 2","Binary Maternal Depression Year 2")
out_var <- c("Telomere length Z-score Year 1","Telomere length Z-score Year 2","Change in Telomere length Z-Score")

tbl3 <- growth_tbl("Maternal Depression", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl3flex <- growth_tbl_flex("Maternal Depression", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl3supp <- growth_tbl("Maternal Depression", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl3flexsupp <- growth_tbl_flex("Maternal Depression", expo_var, out_var, exposure, outcome, H2, H2adj)


#### Table 4 ####
exposure <- c("pss_sum_mom_t3", "pss_sum_dad_t3")
outcome <- c("TS_t3_Z")
expo_var <- c("Maternal Perceived Stress", "Paternal Perceived Stress")
out_var <- c("Telomere Length Z-score Year 2")

tbl4 <- growth_tbl("Parental Stress", expo_var, out_var, exposure, outcome, H3, H3adj, T)
tbl4flex <- growth_tbl_flex("Parental Stress", expo_var, out_var, exposure, outcome, H3, H3adj, T)
tbl4supp <- growth_tbl("Parental Stress", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl4flexsupp <- growth_tbl_flex("Parental Stress", expo_var, out_var, exposure, outcome, H3, H3adj)


#### SAVE TABLES ####

write.csv(tbl2, here('tables/ipv-dep-stress-telo-table1.csv'))
write.csv(tbl3, here('tables/ipv-dep-stress-telo-table2.csv'))
write.csv(tbl4, here('tables/ipv-dep-stress-telo-table3.csv'))

write.csv(tbl2supp, here('tables/ipv-dep-stress-telo-table1supp.csv'))
write.csv(tbl3supp, here('tables/ipv-dep-stress-telo-table2supp.csv'))
write.csv(tbl4supp, here('tables/ipv-dep-stress-telo-table3supp.csv'))

save_as_docx( "Table 1: Association Between Maternal Exposure to IPV and Child Telomere Length" = tbl2flex, 
              "Table 2: Association Between Maternal Depression and Child Telomere Length" = tbl3flex, 
              "Table 3: Association Between Parental Stress and Child Telomere Length" = tbl4flex, 
              path='C:/Users/Sophia/Documents/WASH/WASH IPV and Telomeres/ipv-dep-stress-telo-tables.docx')

save_as_docx( "Table 1: Association Between Maternal Exposure to IPV and Child Telomere Length" = tbl2flexsupp, 
              "Table 2: Association Between Maternal Depression and Child Telomere Length" = tbl3flexsupp, 
              "Table 3: Association Between Parental Stress and Child Telomere Length" = tbl4flexsupp, 
              path='C:/Users/Sophia/Documents/WASH/WASH IPV and Telomeres/ipv-dep-stress-telo-tables supplementary.docx')
