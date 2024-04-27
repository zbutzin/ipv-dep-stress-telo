rm(list=ls())

source(here::here("0-config.R"))
source(here::here("table-functions.R"))

emm <- readRDS(here('results/adjusted/emm_sex_adj_res.RDS'))

tbl1 <- subgroup_tbl("Parental Stress", 
                     c("life_viol_any_t3","viol_any_preg", "viol_any_t2","viol_12m_any_t3", "cesd_sum_t2", "cesd_sum_t2_binary", "cesd_sum_ee_t3", "cesd_sum_ee_t3_binary", "pss_sum_mom_t3", "pss_sum_dad_t3"), 
                     c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z"),
                     c("sex"), 
                     c("Lifetime Exposure to IPV Year 2","Exposure to IPV during Pregnancy", "Exposure to IPV during first year of child's life", "Exposure to IPV in Past 12 Months Year 2", "Maternal Depression Year 1","Binary Maternal Depression Year 1", "Maternal Depression Year 2","Binary Maternal Depression Year 2", "Maternal Perceived Stress", "Paternal Perceived Stress"),
                     c("Telomere length Z-score Year 1", "Telomere length Z-score Year 2","Change in Telomere length Z-Score"),
                     c("Child Sex"), emm)

tbl1

save_as_docx("EMM Table: Effect modification of maternal IPV and child telomere length" = tbl1,
             path = "/Users/kjung0909/Documents/Research/WASHB/ipv-TL/ipv-dep-stress-telo/tables/ipv-telo-emm1.docx",
             pr_section = sect_properties)
