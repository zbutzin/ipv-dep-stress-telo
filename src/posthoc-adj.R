rm(list=ls())

source(here::here("0-config.R"))

d_clean <- readRDS("/Users/kjung0909/Documents/Research/WASHB/bangladesh-cleaned-master-data.RDS") %>% filter(ipv_telo==1)

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth_scaled", "tr")

Wvars[!(Wvars %in% colnames(d_clean))]

#Add in time varying covariates:
#hypothesis Z, outcome=telomere length at year 2
Wvars.HZ.t3<-c(Wvars, "ageday_ht3", "month_ht3", "mhle_month_t3") 
#hypothesis Z, outcome=change in telomere length
Wvars.HZ.d23<-c(Wvars, "ageday_ht2", "ageday_ht3", "month_ht2", "month_ht3", "mhle_month_t3") 
#hypothesis Z, outcome=telomere length at year 1
Wvars.HZ.t2<-c(Wvars, "ageday_ht2", "month_ht2", "mhle_month_t3")

#Loop over exposure-outcome pairs

pick_covariates_HZ <- function(j){
  if(grepl("_t2", j)){Wset = Wvars.HZ.t2}
  if(grepl("_t3", j)){Wset = Wvars.HZ.t3}
  if(grepl("delta", j)){Wset = Wvars.HZ.d23}
  return(Wset)
}

#### Hypothesis Z ####
# Maternal lifetime exposure to type of IPV measured at Year 2 and child telomere length at Year 1, Year 2, and change between Year 1 and Year 2
Xvars <- c("life_emot_viol_t3", "life_phys_viol_t3", "life_sex_viol_t3")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

#Fit models
HZ_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d_clean, X=i, Y=j,  W=pick_covariates_HZ(j))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    HZ_adj_models <- bind_rows(HZ_adj_models, res)
  }
}

#Get primary contrasts
HZ_adj_res <- NULL
for(i in 1:nrow(HZ_adj_models)){
  res <- data.frame(X=HZ_adj_models$X[i], Y=HZ_adj_models$Y[i])
  preds <- predict_gam_diff(fit=HZ_adj_models$fit[i][[1]], d=HZ_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  HZ_adj_res <-  bind_rows(HZ_adj_res , preds$res)
}

#Make list of plots
HZ_adj_plot_list <- NULL
HZ_adj_plot_data <- NULL
for(i in 1:nrow(HZ_adj_models)){
  res <- data.frame(X=HZ_adj_models$X[i], Y=HZ_adj_models$Y[i])
  simul_plot <- gam_simul_CI(HZ_adj_models$fit[i][[1]], HZ_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  HZ_adj_plot_list[[i]] <-  simul_plot$p
  HZ_adj_plot_data <-  rbind(HZ_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(HZ_adj_models, here("models/HZ_adj_models.RDS"))

#Save results
saveRDS(HZ_adj_res, here("results/adjusted/HZ_adj_res.RDS"))

#Save plot data
saveRDS(HZ_adj_plot_data, here("figure-data/HZ_adj_spline_data.RDS"))
