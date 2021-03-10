rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-dep-stress-telo-covariates.RDS"))

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth", "tr")

Wvars[!(Wvars %in% colnames(d))]

#Add in time varying covariates:
#hypothesis 1, outcome=telomere length at year 2
Wvars.H1.t3<-c(Wvars, "ageday_ht3", "month_ht3", "mhle_month_t3") 
#hypothesis 1, outcome=change in telomere length
Wvars.H1.d23<-c(Wvars, "ageday_ht2", "ageday_ht3", "month_ht2", "month_ht3", "mhle_month_t3") 
#hypothesis 1, outcome=telomere length at year 1
Wvars.H1.t2<-c(Wvars, "ageday_ht2", "month_ht2", "mhle_month_t3")
#hypothesis 2, exposure=cesd at year 1, outcome=telomere length at year 1
Wvars2ai<-c(Wvars, "ageday_ht2", "month_ht2", "cesd_month_t2") 
#hypothesis 2, exposure=cesd at year 2, outcome=telomere length at year 2
Wvars2aii<-Wvars.H1.t3
#hypothesis 2, exposure=cesd at year 1, outcome=change in telomere length
Wvars2b<-c(Wvars, "ageday_ht2", "ageday_ht3", "month_ht2", "month_ht3", "cesd_month_t2") 
#hypothesis 2, exposure=cesd at year 1, outcome=telomere length at year 2
Wvars2c<-c(Wvars, "ageday_ht3", "cesd_month_t2", "month_ht3") 
Wvars3mother<-Wvars.H1.t3
Wvars3father<-c(Wvars, "ageday_ht3", "month_ht3", "pss_dad_month_t3")


#Loop over exposure-outcome pairs

pick_covariates_H1 <- function(j){
  if(grepl("_t2", j)){Wset = Wvars.H1.t2}
  if(grepl("_t3", j)){Wset = Wvars.H1.t3}
  if(grepl("delta", j)){Wset = Wvars.H1.d23}
  return(Wset)
}

#### Hypothesis 1 ####
# Maternal exposure to IPV measured at Year 2 is negatively associated with child telomere length
Xvars <- c("life_viol_any_t3")            
Yvars <- c("TS_t3_Z") 

#Fit models
H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates_H1(j))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

Xvars <- c("viol_12m_any_t3")            
Yvars <- c("delta_TS_Z", "TS_t3_Z")

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates_H1(j)) 
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

Xvars <- c("viol_any_preg", "viol_any_t2")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=pick_covariates_H1(j)) 
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  H1_adj_res <-  bind_rows(H1_adj_res , preds$res)
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


#Save models
saveRDS(H1_adj_models, here("models/H1_adj_models.RDS"))

#Save results
saveRDS(H1_adj_res, here("results/adjusted/H1_adj_res.RDS"))

#Save plot data
saveRDS(H1_adj_plot_data, here("figure-data/H1_adj_spline_data.RDS"))




#### Hypothesis 2 ####
Xvars <- c("cesd_sum_t2", "cesd_sum_t2_binary")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")  

#Fit models
H2_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    if(j=="TS_t2_Z"){
      res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars2ai) 
    }else if(j=="delta_TS_Z"){
      res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars2b) 
    }else{
      res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars2c) 
    }
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

Xvars <- c("cesd_sum_ee_t3", "cesd_sum_ee_t3_binary")            
Yvars <- c("TS_t3_Z")  
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars2aii)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

#Get primary contrasts
H2_adj_res <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  if(grepl("binary", H2_adj_models$X[i])){
    preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binaryX=T)
  }else{
    preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }  
  H2_adj_res <-  bind_rows(H2_adj_res , preds$res)
}

#Make list of plots
H2_adj_plot_list <- NULL
H2_adj_plot_data <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_models$fit[i][[1]], H2_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_adj_plot_list[[i]] <-  simul_plot$p
  H2_adj_plot_data <-  rbind(H2_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H2_adj_models, here("models/H2_adj_models.RDS"))

#Save results
saveRDS(H2_adj_res, here("results/adjusted/H2_adj_res.RDS"))

#Save plot data
saveRDS(H2_adj_plot_data, here("figure-data/H2_adj_spline_data.RDS"))




#### Hypothesis 3 ####
#Parental stress measured at Year 2 is negatively associated with child telomere length measured at Year 2. 
Xvars <- c("pss_sum_mom_t3", "pss_sum_dad_t3")            
Yvars <- c("TS_t3_Z") 

#Fit models
H3_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    if(i=="pss_sum_mom_t3"){
      res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars3mother) 
    }else{
      res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wvars3father) 
    }
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_adj_res <-  bind_rows(H3_adj_res , preds$res)
}

#Make list of plots
H3_adj_plot_list <- NULL
H3_adj_plot_data <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H3_adj_models$fit[i][[1]], H3_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_adj_plot_list[[i]] <-  simul_plot$p
  H3_adj_plot_data <-  rbind(H3_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
saveRDS(H3_adj_models, here("models/H3_adj_models.RDS"))

#Save results
saveRDS(H3_adj_res, here("results/adjusted/H3_adj_res.RDS"))

#Save plot data
saveRDS(H3_adj_plot_data, here("figure-data/H3_adj_spline_data.RDS"))

