rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(patchwork)

#load spline data
H1_spline <- readRDS(here("figure-data/H1_adj_spline_data.RDS"))

#load results for quartiles
H1_quartiles <- readRDS(here("results/adjusted/H1_adj_res.RDS"))

d_for_plot <- function(x_name, y_name, x_var, y_var, spline, quart){
  d <- NULL
  for (i in 1: length(x_var)) {
    for (j in 1:length(y_var)){
      exists <- (quart%>%filter(X==x_var[i], Y==y_var[j]) %>% nrow()) != 0
      if (exists){
        new <- data.frame(x=x_name[i], y=y_name[j], quart%>%filter(X==x_var[i], Y==y_var[j]))
        d <- rbind(d, new)
      }
    }
  }
  d
}


d1 <- d_for_plot(c("During lifetime","During Pregnancy", "During first year of child's life", "During second year of child's life"),
                 c("TL Z-score Year 1", "TL Z-score Year 2","Change in TL Z-Score"), 
                 c("life_viol_any_t3","viol_any_preg", "viol_any_t2","viol_12m_any_t3"),   
                 c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z"),
                 H1_spline, H1_quartiles) #%>% left_join(colors,"x")

d1$x <- factor(d1$x,levels=c("During lifetime","During Pregnancy", "During first year of child's life", "During second year of child's life"))
#d1$x <- factor(d1$x, levels=rev(levels(d1$x)))
d1$y <- factor(d1$y)
d1$`Time of outcome measurement` <- factor(ifelse(grepl("t3", d1$Y), "Year 2", "Year 1"))

# Forest Plot
p <- ggplot(d1, aes(x=x, y=point.diff)) + 
  geom_pointrange(aes(ymin=lb.diff , ymax=ub.diff, color=X),
                  position = position_dodge2(width = 0.5),
                  size = 1, show.legend = F) +
  facet_grid(rows=vars(y), scales = "free") + 
  labs(y = "Adjusted difference in mean child telomere length outcome\nbetween absence and presence of maternal exposure", 
       x="Maternal exposure to intimate partner violence") +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_shape_manual(breaks=c("Year 1","Year 2"), values=c(21,16)) +
  guides(color="none")+
  scale_color_manual(values = tableau10[c(7:1)]) + 
  theme_ki() +
  theme(plot.title = element_text(hjust = 0),
        axis.title.x = element_text(size=11),
        strip.text = element_text(hjust=0.5, size=10),
        strip.placement = "outside",
        axis.text.y = element_text(hjust = 1, size=9),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "bottom")      
p


p %>% ggsave(filename="figures/ipv-telo_point_diff_ipv.jpg", width=10, height=7)

# Spline plots
# not useful since exposure is binary
d_for_plot <- function(x_name, y_name, x_var, y_var, spline, quart){
  d <- NULL
  for (i in 1: length(x_var)) {
    for (j in 1:length(y_var)){
      exists <- (quart%>%filter(X==x_var[i], Y==y_var[j]) %>% nrow()) != 0
      if (exists){
        new <- data.frame(x=x_name[i], y=y_name[j], spline%>%filter(Xvar==x_var[i], Yvar==y_var[j]), quart%>%filter(X==x_var[i], Y==y_var[j])%>%select(q1, q3))
        d <- rbind(d, new)
      }
    }
  }
  d
}

d1 <- d_for_plot(c("Lifetime Exposure to IPV Year 2","Exposure to IPV during Pregnancy",
                   "Exposure to IPV during first year of child's life", "Exposure to IPV in Past 12 Months Year 2"),
                 c("Telomere length Z-score Year 1", "Telomere length Z-score Year 2","Change in Telomere length Z-Score"), 
                 c("life_viol_any_t3","viol_any_preg", "viol_any_t2","viol_12m_any_t3"),   
                 c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z"),
                 H1_spline, H1_quartiles) #%>% left_join(colors,"x")

d1$x <- factor(d1$x,levels=c("Lifetime Exposure to IPV Year 2","Exposure to IPV during Pregnancy",
                             "Exposure to IPV during first year of child's life", "Exposure to IPV in Past 12 Months Year 2"))
#d1$x <- factor(d1$x, levels=rev(levels(d1$x)))
d1$y <- factor(d1$y)

d1 <- d1 %>% filter(!grepl("def", Xvar))

t2 <- d1 %>% filter(grepl("t2", Yvar)) %>% ggplot(aes(x=X))+
  geom_smooth(aes(y = fit, color=x), se = F) +
  geom_vline(aes(xintercept=q1), size=.5, color="grey30", linetype="dashed") +
  geom_vline(aes(xintercept=q3), size=.5, color="grey30", linetype="dashed") +
  #geom_point(aes(y=Y), alpha=0.5) +
  geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=x, color=x), alpha=0.5) + 
  scale_colour_manual(values=tableau10[c(1,5,6,4)]) + 
  scale_fill_manual(values=tableau10[c(1,5,6,4)]) + 
  xlab(" ") + ylab(element_blank()) + 
  facet_grid(cols = vars(x), rows = vars(y), scales = "free") +
  theme_ki() +
  theme(strip.text.x = element_text(size=10),
        strip.text.y = element_text(size=10),
        panel.spacing = unit(.3, "lines"))      
t2

t3 <- d1 %>% filter(grepl("t3", Yvar)) %>% ggplot(aes(x=X))+
  geom_smooth(aes(y = fit, color=x), se = F) +
  geom_vline(aes(xintercept=q1), size=.5, color="grey30", linetype="dashed") +
  geom_vline(aes(xintercept=q3), size=.5, color="grey30", linetype="dashed") +
  #geom_point(aes(y=Y), alpha=0.5) +
  geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=x, color=x), alpha=0.5) + 
  scale_colour_manual(values=tableau10[c(1,5,6,4)]) + 
  scale_fill_manual(values=tableau10[c(1,5,6,4)]) + 
  xlab(" ") + ylab(element_blank()) + 
  facet_grid(cols = vars(x), rows = vars(y), scales="free") +
  theme_ki() +
  theme(strip.text.x = element_text(size=10),
        strip.text.y = element_text(hjust=0.5, size=10),
        panel.spacing = unit(.3, "lines"))      

t3
H1_quartiles

t2 %>% ggsave(filename="figures/adj-splines-t2.jpg", width = 8, height = 9)
t3 %>% ggsave(filename="figures/adj-splines-t3.jpg", width = 8, height = 9)
