rm(list=ls())

library('flextable')
library('officer')
source(here::here("0-config.R"))
library('here')
library('data.table')

#install.packages("boxr")
#library(boxr)
#usethis::edit_r_environ()
#d <-box_read("871638120165") # %>% filter(immune_dev==1)
#d <- readRDS("/Users/lgg/Box/washb/Bangladesh/Master\ Dataset/bangladesh-cleaned-master-data.RDS") %>% filter(ipv_telo==1)
d_clean <- readRDS("/Users/kjung0909/Documents/Research/WASHB/bangladesh-cleaned-master-data.RDS") %>% filter(ipv_telo==1)

# load enrollment characteristics and results
HZ <- readRDS(here('results/unadjusted/HZ_res.RDS'))
HZadj <- readRDS(here('results/adjusted/HZ_adj_res.RDS'))

#### Functions for growth tables ####
source(here::here("table-functions.R"))


#### MAIN TABLES ####
#### Table Z ####
# type of IPV

exposure <- c("life_emot_viol_t3", "life_phys_viol_t3", "life_sex_viol_t3")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Lifetime Exposure to Emotional Violence Year 2","Lifetime Exposure to Physical Violence Year 2",
              "Lifetime Exposure to Sexual Violence Year 2")
out_var <- c("Telomere length Z-score Year 1", "Telomere length Z-score Year 2","Change in Telomere length Z-Score")

tblZ <- growth_tbl("Type of IPV", expo_var, out_var, exposure, outcome, HZ, HZadj, T)
tblZflex <- growth_tbl_flex("Type of IPV", expo_var, out_var, exposure, outcome, HZ, HZadj, T, 1.6, 1.3)
#tblZsupp <- growth_tbl("Type of IPV", expo_var, out_var, exposure, outcome, HZ, HZadj)
#tblZflexsupp <- growth_tbl_flex("Type of IPV", expo_var, out_var, exposure, outcome, HZ, HZadj)


#### SAVE TABLES ####

write.csv(tblZ, here('tables/ipv-dep-stress-telo-tableZ.csv'))

#write.csv(tblZsupp, here('tables/ipv-dep-stress-telo-tableZsupp.csv'))

save_as_docx( "Table 5: Association Between Maternal Exposure to Type of IPV and Child Telomere Length" = tblZflex, 
              path='/Users/kjung0909/Documents/Research/WASHB/ipv-TL/ipv-dep-stress-telo/tables/ipv-telo-posthoc-tables.docx',
              pr_section = sect_properties)

#save_as_docx( "Table 5: Association Between Maternal Exposure to Type of IPV and Child Telomere Length" = tblZflexsupp, path='/Users/kjung0909/Documents/Research/WASHB/ipv-TL/ipv-dep-stress-telo/tables/ipv-dep-stress-telo-tables supplementary.docx')


#### Count exposure by type of IPV ####
# Unique mothers
exp <- c("viol_any_t2", "life_viol_any_t3", "cesd_sum_t2", "life_viol_any_t3","viol_12m_any_t3", "pss_sum_mom_t3", "pss_sum_dad_t3", "cesd_sum_ee_t3") 
out <- c("TS_t2","TS_t3", "delta_TS") 

filtering <- function(row){
  any(!is.na(row))
}

m <- d_clean[apply(select(d_clean, all_of(out)), 1, filtering),]

sum(m$life_phys_viol_t3, na.rm=T)
sum(m$life_emot_viol_t3, na.rm=T)
sum(m$life_sex_viol_t3, na.rm=T)
