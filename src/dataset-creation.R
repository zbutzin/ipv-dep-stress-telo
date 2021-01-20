rm(list=ls())
source(here::here("0-config.R"))

#load maternal exposure data
#load(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.RData"))

# load telo-growth data (also includes maternal exposure)
d<-read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-telo-growth-covariates-telolab-anthro.csv"))
names(d)

#Z-scores of telomere measurements
ipv.dep.stress.telo <- d %>% 
   mutate(TS_t2_Z = scale(TS_t2, center=TRUE, scale=TRUE)[,1]) %>%
   mutate(TS_t3_Z = scale(TS_t3, center=TRUE, scale=TRUE)[,1]) %>%
   mutate(delta_TS_Z = scale(delta_TS, center=TRUE, scale=TRUE)[,1])

#recode mothers who responded no to lifetime ipv as 0 for ipv in the last 12 mo
ipv.dep.stress.telo$viol_12m_any_t3_recode <- ifelse(is.na(ipv.dep.stress.telo$viol_12m_any_t3) & 
                                                       ipv.dep.stress.telo$life_viol_any_t3 == 0, 
                                                     0, ipv.dep.stress.telo$viol_12m_any_t3)

# add variables to turn cesd into binary variables
# classify top 25% of mothers in sample as experiencing high depressive symptoms
cesd_t2_q<-quantile(ipv.dep.stress.telo$cesd_sum_t2, na.rm=T)[4]
cesd_t3_q<-quantile(ipv.dep.stress.telo$cesd_sum_ee_t3, na.rm=T)[4]
ipv.dep.stress.telo$cesd_sum_t2_binary<-ifelse(ipv.dep.stress.telo$cesd_sum_t2 >= cesd_t2_q, 
                                               1, 0)
ipv.dep.stress.telo$cesd_sum_ee_t3_binary<-ifelse(ipv.dep.stress.telo$cesd_sum_ee_t3 >= cesd_t3_q, 
                                               1, 0)

# add hh wealth 
hhwealth <- read.csv("C:/Users/Sophia/Documents/ee-secondary/sophia scripts/hhwealth.csv")
full <- left_join(ipv.dep.stress.telo, hhwealth, by="dataid")

saveRDS(full, paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-dep-stress-telo-covariates.RDS"))
