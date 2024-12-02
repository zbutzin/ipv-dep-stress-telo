#Correlation of TS2 and TS3 and reviewer response

rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/ipv-dep-stress-telo-covariates.RDS"))

t2t3cor <- cor.test(d$TS_t3_Z, d$TS_t2_Z)

#correlations for reviewer response:
# What was the correlation between maternal IPV, depression and parental PSS?
# Are their correlations consistent with the effect direction with telomere length?
print_cor <- function(x,y, xlab, ylab){
  return(data.frame(X=xlab, Y=ylab, r=round(cor(d[[x]], d[[y]], use = "pairwise.complete.obs"),3)))
}

cor_tab <- bind_rows(
  print_cor("life_viol_any_t3", "cesd_sum_ee_t3", "Lifetime maternal IPV", "Maternal depression (CESD)"),
  print_cor("life_viol_any_t3", "pss_sum_mom_t3", "Lifetime maternal IPV", "Maternal stress  (PSS)"),
  print_cor("life_viol_any_t3", "pss_sum_dad_t3", "Lifetime maternal IPV", "Paternal stress (PSS)"),
  print_cor("cesd_sum_ee_t3", "pss_sum_mom_t3", "Maternal depression (CESD)", "Maternal stress  (PSS)"),
  print_cor("cesd_sum_ee_t3", "pss_sum_dad_t3", "Maternal depression (CESD)", "Paternal stress (PSS)"),
  print_cor("life_viol_any_t3", "delta_TS_Z", "Lifetime maternal IPV", "Change in Telomere Length"),
  print_cor("cesd_sum_ee_t3", "delta_TS_Z", "Maternal depression (CESD)", "Change in Telomere Length"),
  print_cor("pss_sum_mom_t3", "delta_TS_Z", "Maternal stress  (PSS)", "Change in Telomere Length"),
  print_cor("pss_sum_dad_t3", "delta_TS_Z", "Paternal stress (PSS)", "Change in Telomere Length"),
  print_cor("life_viol_any_t3", "TS_t2_Z", "Lifetime maternal IPV", "Telomere Length - T2"),
  print_cor("cesd_sum_ee_t3", "TS_t2_Z", "Maternal depression (CESD)", "Telomere Length - T2"),
  print_cor("pss_sum_mom_t3", "TS_t2_Z", "Maternal stress  (PSS)", "Telomere Length - T2"),
  print_cor("pss_sum_dad_t3", "TS_t2_Z", "Paternal stress (PSS)", "Telomere Length - T2"),
  print_cor("life_viol_any_t3", "TS_t3_Z", "Lifetime maternal IPV", "Telomere Length - T3"),
  print_cor("cesd_sum_ee_t3", "TS_t3_Z", "Maternal depression (CESD)", "Telomere Length - T3"),
  print_cor("pss_sum_mom_t3", "TS_t3_Z", "Maternal stress  (PSS)", "Telomere Length - T3"),
  print_cor("pss_sum_dad_t3", "TS_t3_Z", "Paternal stress (PSS)", "Telomere Length - T3")
)

knitr::kable(cor_tab)
