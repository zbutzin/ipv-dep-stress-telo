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


# check covariate missingness and categorize momheight
Wvars<-c("sex", "birthord", "momage","momheight","momedu", 
         "hfiacat", "Ncomp", "watmin", "walls", "floor", "tr", 
         "HHwealth", "ageday_ht2", "ageday_ht3", "month_ht2", "cesd_month_t2", 
         "month_ht3", "mhle_month_t3", "pss_dad_month_t3")

W <- full %>% select(all_of(Wvars))  

miss <- data.frame(name = names(W), missing = colSums(is.na(W)), row.names = c(1:ncol(W)))
for (i in 1:nrow(miss)) {
  miss$class[i] <- class(W[,which(colnames(W) == miss[i, 1])])
}
miss

for(outcome in c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")){
  d_sub <- subset(full, !is.na(full[,outcome]))
  W_sub <- d_sub %>% select(all_of(Wvars))  
  
  miss_sub <- data.frame(name = names(W_sub), missing = colSums(is.na(W_sub)), row.names = c(1:ncol(W_sub)))
  for (i in 1:nrow(miss_sub)) {
    miss_sub$class[i] <- class(W_sub[,which(colnames(W_sub) == miss_sub[i, 1])])
  }
  print(outcome)
  print(miss_sub)
}

saveRDS(full, paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-dep-stress-telo-covariates.RDS"))
