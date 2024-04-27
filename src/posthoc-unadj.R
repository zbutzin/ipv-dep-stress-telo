rm(list=ls())

library('haven')
library('flextable')
library('officer')
source(here::here("0-config.R"))
library('here')
library('data.table')

#d_master <- read_dta('MHLE_Mother_clean_data_30May16_deidentified2.dta')
d_clean <- readRDS("/Users/kjung0909/Documents/Research/WASHB/bangladesh-cleaned-master-data.RDS") %>% filter(ipv_telo==1)

# checking NA in type of ipv
#sum(is.na(d_clean$emot_viol_preg))
#length(d_clean$emot_viol_preg)

#col_names <- c('emot_viol_preg', 'sex_viol_preg', 'phys_viol_preg', 'emot_viol_t2', 'sex_viol_t2', 'phys_viol_t2', 'emot_viol_12m_t3', 'sex_viol_12m_t3', 'phys_viol_12m_t3', 'emot_viol_t3', 'sex_viol_t3', 'phys_viol_t3')

#for (i in col_names) {
#    column_vector_name <- paste0("d_clean$", i)
#    column_vector <- d_clean[, i]
#      length <- length(column_vector)
#      num_na <- sum(is.na(column_vector))
#      percentage <- round(num_na/length*100, 1)
#      print(paste0("variable: ", i, "; NAs: ", num_na, "; total values: ", length, "; percentage: ", percentage))
#}

#rm(list=ls())

#d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-dep-stress-telo-covariates.RDS"))
#names(d)

#Example:

#Fit GAM model with random effects for childid
#res_unadj <- fit_RE_gam(d=d, X="t3_cort_z01", Y="laz_t3",  W=NULL)

#Get predictions of differences from the 25th percentile of exposure
#preds_unadj <- predict_gam_diff(fit=res_unadj$fit, d=res_unadj$dat, quantile_diff=c(0.25,0.75), Xvar="delta_TS", Yvar="laz_t3")


#Primary parameter we are estimating: difference between 25th and 75th percentile of the exposure
#preds_unadj$res

#Plot the difference from the 25th percentile for the full range of the exposure:
#NOTE: not making these plots anymore, just using for diagnostics
#p <- plot_gam_diff(preds_unadj$plotdf)
#print(p)

#Fit spline with simultaneous confidence intervals
#simul_plot <- gam_simul_CI(res_unadj$fit, res_unadj$dat, xlab="delta_TS", ylab="laz_t3", title="example title")
#simul_plot$p


#Loop over exposure-outcome pairs

#### Hypothesis Z: Post-hoc type of IPV ####
# Maternal exposure to cumulative lifetime type of IPV measured at Year 2 with child telomere length measured at Year 1, Year 2, and change between Year 1 and Year 2
Xvars <- c("life_emot_viol_t3", "life_phys_viol_t3", "life_sex_viol_t3")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

#Fit models
HZ_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d_clean, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    HZ_models <- bind_rows(HZ_models, res)
  }
}

#Get primary contrasts
HZ_res <- NULL
for(i in 1:nrow(HZ_models)){
  res <- data.frame(X=HZ_models$X[i], Y=HZ_models$Y[i])
  preds <- predict_gam_diff(fit=HZ_models$fit[i][[1]], d=HZ_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  HZ_res <-  bind_rows(HZ_res , preds$res)
}

#Make list of plots
HZ_plot_list <- NULL
HZ_plot_data <- NULL
for(i in 1:nrow(HZ_models)){
  res <- data.frame(X=HZ_models$X[i], Y=HZ_models$Y[i])
  simul_plot <- gam_simul_CI(HZ_models$fit[i][[1]], HZ_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  HZ_plot_list[[i]] <-  simul_plot$p
  HZ_plot_data <-  rbind(HZ_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


saveRDS(HZ_res, file=here('results/unadjusted/HZ_res.RDS'))
saveRDS(HZ_models, file=here('models/HZ_models.RDS'))
saveRDS(HZ_plot_data, file=here('figure-data/HZ_unadj_spline_data.RDS'))




