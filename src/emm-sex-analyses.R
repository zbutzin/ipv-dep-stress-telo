rm(list=ls())

source(here::here("0-config.R"))
source(here::here("table-functions.R"))

d <- readRDS("/Users/kjung0909/Documents/Research/WASHB/bangladesh-cleaned-master-data.RDS") %>% filter(.$ipv_telo == 1)
#d<-readRDS(paste0(dropboxDir, "Data/Cleaned/Audrie/pregnancy_child_immune_covariates_data.RDS"))

#Maternal ipv/dep/stress and child telo EMM by child sex
Xvars <- c("life_viol_any_t3","viol_any_preg", "viol_any_t2","viol_12m_any_t3", "cesd_sum_t2", "cesd_sum_t2_binary", "cesd_sum_ee_t3", "cesd_sum_ee_t3_binary", "pss_sum_mom_t3", "pss_sum_dad_t3")            
Yvars <- c("TS_t2", "TS_t3", "delta_TS")

#Fit models
H1_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL, V="sex")
    res <- data.frame(X=i, Y=j, V="sex", int.p =res_unadj$int.p, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}

#Debug the gam_EMM
H1_models$X

#Get primary contrasts
H1_res <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  if(grepl("viol|binary", H1_models$X[i])){
    preds <- predict_gam_emm(fit=H1_models$fit[i][[1]], d=H1_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=H1_models$X[i], Yvar=H1_models$Y[i], binaryX=T)
  }else{
    preds <- predict_gam_emm(fit=H1_models$fit[i][[1]], d=H1_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=H1_models$X[i], Yvar=H1_models$Y[i])
  }
  gamm_diff_res <- data.frame(V=H1_models$V[i] , preds$res) %>% mutate(int.Pval = c(NA, H1_models$int.p[[i]]))
  
  H1_res <-  bind_rows(H1_res , gamm_diff_res)
}

#Make list of plots
H1_plot_list <- NULL
H1_plot_data <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  simul_plot <- gam_simul_CI(H1_models$fit[i][[1]], H1_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_plot_list[[i]] <-  simul_plot$p
  H1_plot_data <-  rbind(H1_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}

H1_res <- H1_res %>% mutate(BH.Pval=p.adjust(Pval, method="BH"),
                            BH.Pval.int=p.adjust(int.Pval, method="BH")) 

#Save models
saveRDS(H1_models, here("models/emm_sex.RDS"))

#Save results
saveRDS(H1_res, here("results/unadjusted/emm_sex_res.RDS"))

#Save plot data
#saveRDS(H1_plot_data, here("figure-data/cytokine-ratios_unadj_spline_data.RDS"))

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("birthord", "momage","momheight","momedu","gest_age_weeks", "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth_scaled")

Wvars[!(Wvars %in% colnames(d))]

#Add in time varying covariates:
Wvars2 <- c(Wvars, c("ageday_bt2", "month_blood_t0", "month_bt2"))
Wvars3 <- c(Wvars, c("ageday_bt3", "month_blood_t0", "month_bt3"))

pick_covariates <- function(j){
  # j is outcome as string
  # choose correct adjustment set based on outcome
  if(grepl("t2", j)){Wset = Wvars2}
  else{Wset = Wvars3}
  return(Wset)
}

H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset, V="sex")
    res <- data.frame(X=i, Y=j, V="sex", int.p =res_adj$int.p, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  if(grepl("_def", H1_adj_models$X[i])){
    preds <- predict_gam_emm(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=H1_adj_models$X[i], Yvar=H1_adj_models$Y[i], binaryX=T)
  }else{
    preds <- predict_gam_emm(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=H1_adj_models$X[i], Yvar=H1_adj_models$Y[i])
  }
  gamm_diff_res <- data.frame(V=H1_adj_models$V[i] , preds$res) %>% mutate(int.Pval = c(NA, H1_adj_models$int.p[[i]]))
  
#  H1_adj_res <-  bind_rows(H1_adj_res , gamm_diff_res)
}


#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}

H1_adj_res <- H1_adj_res %>% mutate(BH.Pval=p.adjust(Pval, method="BH"),
                                    BH.Pval.int=p.adjust(int.Pval, method="BH")) 

#Save results
saveRDS(H1_adj_res, here("results/adjusted/emm_sex_adj_res.RDS"))

#Save plot data
saveRDS(H1_adj_plot_data, here("figure-data/emm_sex_adj_spline.data.RDS"))

