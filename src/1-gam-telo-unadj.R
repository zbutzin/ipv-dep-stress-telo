rm(list=ls())

source(here::here("0-config.R"))

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

#### Hypothesis 1a ####
# Maternal exposure to cumulative lifetime IPV measured at Year 2 is negatively associated with child telomere length measured at Year 2
Xvars <- c("life_viol_any_t3")            
Yvars <- c("TS_t3_Z") 

#Fit models
H1_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}

#### Hypothesis 1b-c ####
# Maternal exposure to IPV from Year 1 to Year 2 is positively associated with child telomere shortening from Year 1 to Year 2
# Maternal exposure to IPV from Year 1 to Year 2 is negatively associated with child telomere length at Year 2
Xvars <- c("viol_12m_any_t3")            
Yvars <- c("delta_TS_Z", "TS_t3_Z") 

#Fit models
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}


#### Hypothesis 1d-i ####
# Maternal exposure to IPV during pregnancy is negatively associated with child telomere length
# Maternal exposure to IPV from Year 0 to Year 1 is negatively associated with child telomere length
Xvars <- c("viol_any_preg", "viol_any_t2")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

#Fit models
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}

#Get primary contrasts
H1_res <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  preds <- predict_gam_diff(fit=H1_models$fit[i][[1]], d=H1_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  H1_res <-  bind_rows(H1_res , preds$res)
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


saveRDS(H1_res, file=here('results/unadjusted/H1_res.RDS'))
saveRDS(H1_models, file=here('models/H1_models.RDS'))
saveRDS(H1_plot_data, file=here('figure-data/H1_unadj_spline_data.RDS'))



#### Hypothesis 2a ####
#Maternal depression measured at Years 1 and 2 is negatively associated with concurrent child telomere length at Years 1 and 2
Xvars <- c("cesd_sum_t2", "cesd_sum_ee_t3")            
Yvars <- c("TS_t2_Z", "TS_t3_Z") 

#Fit models
H2a_models <- NULL
for(i in c(1:2)){
  res_unadj <- fit_RE_gam(d=d, X=Xvars[i], Y=Yvars[i],  W=NULL)
  res <- data.frame(X=Xvars[i], Y=Yvars[i], fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
  H2a_models <- bind_rows(H2a_models, res)
}

Xvars <- c("cesd_sum_t2_binary", "cesd_sum_ee_t3_binary")
Yvars <- c("TS_t2_Z", "TS_t3_Z") 
for(i in c(1:2)){
  res_unadj <- fit_RE_gam(d=d, X=Xvars[i], Y=Yvars[i],  W=NULL)
  res <- data.frame(X=Xvars[i], Y=Yvars[i], fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
  H2a_models <- bind_rows(H2a_models, res)
}

#Get primary contrasts
H2a_res <- NULL
for(i in 1:nrow(H2a_models)){
  res <- data.frame(X=H2a_models$X[i], Y=H2a_models$Y[i])
  if(grepl("binary", H2a_models$X[i])){
    preds <- predict_gam_diff(fit=H2a_models$fit[i][[1]], d=H2a_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  }else{
    preds <- predict_gam_diff(fit=H2a_models$fit[i][[1]], d=H2a_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }
  H2a_res <-  bind_rows(H2a_res , preds$res)
}

#Make list of plots
H2a_plot_list <- NULL
H2a_plot_data <- NULL
for(i in 1:nrow(H2a_models)){
  res <- data.frame(X=H2a_models$X[i], Y=H2a_models$Y[i])
  simul_plot <- gam_simul_CI(H2a_models$fit[i][[1]], H2a_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2a_plot_list[[i]] <-  simul_plot$p
  H2a_plot_data <-  rbind(H2a_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H2a_models, here("models/H2a_models.RDS"))

#Save results
saveRDS(H2a_res, here("results/unadjusted/H2a_res.RDS"))


#Save plots
#saveRDS(H2a_plot_list, here("figure-objects/H2a_unadj_splines.RDS"))

#Save plot data
saveRDS(H2a_plot_data, here("figure-data/H2a_unadj_spline_data.RDS"))


#### Hypothesis 2b ####
#Maternal depression at Year 1 is positively associated with child telomere shortening between Year 1 and Year 2
Xvars <- c("cesd_sum_t2", "cesd_sum_t2_binary")            
Yvars <- c("delta_TS_Z") 

#Fit models
H2b_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H2b_models <- bind_rows(H2b_models, res)
  }
}

#Get primary contrasts
H2b_res <- NULL
for(i in 1:nrow(H2b_models)){
  res <- data.frame(X=H2b_models$X[i], Y=H2b_models$Y[i])
  if(grepl("binary", H2b_models$X[i])){
    preds <- predict_gam_diff(fit=H2b_models$fit[i][[1]], d=H2b_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  }else{
    preds <- predict_gam_diff(fit=H2b_models$fit[i][[1]], d=H2b_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }
  H2b_res <-  bind_rows(H2b_res , preds$res)
}

#Make list of plots
H2b_plot_list <- NULL
H2b_plot_data <- NULL
for(i in 1:nrow(H2b_models)){
  res <- data.frame(X=H2b_models$X[i], Y=H2b_models$Y[i])
  simul_plot <- gam_simul_CI(H2b_models$fit[i][[1]], H2b_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2b_plot_list[[i]] <-  simul_plot$p
  H2b_plot_data <-  rbind(H2b_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H2b_models, here("models/H2b_models.RDS"))

#Save results
saveRDS(H2b_res, here("results/unadjusted/H2b_res.RDS"))


#Save plots
#saveRDS(H2b_plot_list, here("figure-objects/H2b_unadj_splines.RDS"))

#Save plot data
saveRDS(H2b_plot_data, here("figure-data/H2b_unadj_spline_data.RDS"))



#### Hypothesis 2c ####
#Maternal depression at Year 1 is negatively associated with subsequent child telomere length at Year 2
Xvars <- c("cesd_sum_t2", "cesd_sum_t2_binary")            
Yvars <- c("TS_t3_Z") 

#Fit models
H2c_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H2c_models <- bind_rows(H2c_models, res)
  }
}

#Get primary contrasts
H2c_res <- NULL
for(i in 1:nrow(H2c_models)){
  res <- data.frame(X=H2c_models$X[i], Y=H2c_models$Y[i])
  if(grepl("binary", H2c_models$X[i])){
    preds <- predict_gam_diff(fit=H2c_models$fit[i][[1]], d=H2c_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  }else{
    preds <- predict_gam_diff(fit=H2c_models$fit[i][[1]], d=H2c_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }
  H2c_res <-  bind_rows(H2c_res , preds$res)
}

#Make list of plots
H2c_plot_list <- NULL
H2c_plot_data <- NULL
for(i in 1:nrow(H2c_models)){
  res <- data.frame(X=H2c_models$X[i], Y=H2c_models$Y[i])
  simul_plot <- gam_simul_CI(H2c_models$fit[i][[1]], H2c_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2c_plot_list[[i]] <-  simul_plot$p
  H2c_plot_data <-  rbind(H2c_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H2c_models, here("models/H2c_models.RDS"))

#Save results
saveRDS(H2c_res, here("results/unadjusted/H2c_res.RDS"))


#Save plots
#saveRDS(H2c_plot_list, here("figure-objects/H2c_unadj_splines.RDS"))

#Save plot data
saveRDS(H2c_plot_data, here("figure-data/H2c_unadj_spline_data.RDS"))

#### Hypothesis 3 ####
#Parental stress measured at Year 2 is negatively associated with child telomere length measured at Year 2. 
Xvars <- c("pss_sum_mom_t3", "pss_sum_dad_t3")            
Yvars <- c("TS_t3_Z") 

#Fit models
H3_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H3_models <- bind_rows(H3_models, res)
  }
}

#Get primary contrasts
H3_res <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  preds <- predict_gam_diff(fit=H3_models$fit[i][[1]], d=H3_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_res <-  bind_rows(H3_res , preds$res)
}

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  simul_plot <- gam_simul_CI(H3_models$fit[i][[1]], H3_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


#Save models
saveRDS(H3_models, here("models/H3_models.RDS"))

#Save results
saveRDS(H3_res, here("results/unadjusted/H3_res.RDS"))


#Save plots
#saveRDS(H3_plot_list, here("figure-objects/H3_unadj_splines.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("figure-data/H3_unadj_spline_data.RDS"))


# save overall hypothesis results/models/figure data
# H2
H2_res <- rbind(H2a_res, H2b_res, H2c_res)
H2_models <- rbind(H2a_models, H2b_models, H2c_models)
H2_plot_data <-rbind(H2a_plot_data, H2b_plot_data, H2c_plot_data)

saveRDS(H2_res, file=here('results/unadjusted/H2_res.RDS'))
saveRDS(H2_models, file=here('models/H2_models.RDS'))
saveRDS(H2_plot_data, file=here('figure-data/H2_unadj_spline_data.RDS'))


