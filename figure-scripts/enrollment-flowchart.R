rm(list=ls())
source(here::here("0-config.R"))
library(tibble)
data <- tibble(x = -10:100, y= -10:100)
head(data)

d <- readRDS("/Users/sophiatan/Library/CloudStorage/Box-Box/washb/Bangladesh/Master Dataset/bangladesh-cleaned-master-data.RDS") %>% filter(ipv_telo==1)
exposures_y1 <- c("viol_any_t2", "life_viol_any_t3", "cesd_sum_t2")
outcomes_y1 <- c("TS_t2") 
exposures_y2 <- c("life_viol_any_t3","viol_12m_any_t3", "pss_sum_mom_t3", "pss_sum_dad_t3", "cesd_sum_ee_t3")
outcomes_y2 <- c("TS_t3", "delta_TS")

#function for filtering for only participants with at least one outcome
filtering <- function(row){
  any(!is.na(row))}

y1_has_exposures<-d[apply(select(d, all_of(exposures_y1)), 1, filtering),]
y1_has_both<-y1_has_exposures[apply(select(y1_has_exposures, all_of(outcomes_y1)), 1, filtering),]
y1_has_both <- y1_has_both %>% group_by(tr) %>% summarise(clusterid=length(unique(clusterid)),
                                           n=n())


if(is.null(exposures_y2)){
  y2_has_exposures <- y1_has_exposures
}else{
  y2_has_exposures <- d[apply(select(d, all_of(exposures_y2)), 1, filtering),]
}
y2_has_both<-y2_has_exposures[apply(select(y2_has_exposures, all_of(outcomes_y2)), 1, filtering),]
y2_has_both <- y2_has_both %>% group_by(tr) %>% summarise(clusterid=length(unique(clusterid)),
                                                          n=n())


library(tibble)
data <- tibble(x = 1:100, y= 1:100)
head(data)

library(dplyr)
data %>% 
  ggplot(aes(x, y)) +
  scale_x_continuous(minor_breaks = seq(10, 100, 10)) +
  scale_y_continuous(minor_breaks = seq(10, 100, 10)) +
  theme_void() ->
  p

p +
  geom_rect(xmin = 25, xmax=75, ymin=100, ymax=104, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=102,label= '13,279 compounds assessed for eligibility', size=3.5) ->#+
  #annotate('text', x= 50, y=102,label= 'Figure S1: CONSORT Diagram for the WASH Benefits immune status and growth factor study population', size=3) ->
  p

p +
  geom_rect(xmin = 58, xmax=104, ymin=91, ymax=99, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 81, y=95,label= 'Excluded: 7,728 compounds \n 7,429 compounds excluded to create buffer zones\n 219 compounds did not meet enrollment criteria\n 80 compounds declined to participate', size=3.5) +
  annotate('text', x= 3, y=95,label= 'Enrollment', size=4) +
  geom_rect(xmin = 25, xmax=75, ymin=83.5, ymax=89.5, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 50, y=87.5,label= '
720 clusters created and randomly allocated across 7 arms \n 5,551 compounds randomly allocated across 7 arms \n 2 out of 7 arms selected into substudy', size=3.5)  +
  annotate('text', x= 2.5, y=86.5,label= 'Allocation', size=4) +
  geom_rect(xmin = 9, xmax=25, ymin=76, ymax=82, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 17, y=80,label= '
Control \n 180 clusters \n 1,382 households', size=3.5) +
  geom_rect(xmin = 71, xmax=104, ymin=76, ymax=82, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 88, y=80,label= '
Nutrition+Water+Sanitation+Handwashing \n 90 clusters \n 686 households ', size=3.5) +
  geom_rect(xmin = 71, xmax=104, ymin=63, ymax=75, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 88, y=70,label= '
Year 1 \n 63 clusters \n 480 children \n Year 2 \n 67 clusters \n 505 children ', size=3.5)+
  geom_rect(xmin = 71, xmax=104, ymin=32, ymax=62, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 88, y=48,label= '
Year 1 \n 100 children lost to follow-up \n 9 moved \n 29 absent \n 14 withdrew \n 37 no live birth \n 11 child death \n Year 2 \n 25 new children measured  \n 104 children lost to follow-up \n 28 moved \n 2 absent \n 18 withdrew \n 38 no live birth \n 18 child death ', size=3.5) +
  geom_rect(xmin = 71, xmax=104, ymin=19, ymax=31, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 88, y=26,label= '
Year 1 \n 63 clusters \n 380 children \n Year 2 \n 67 clusters \n 401 children ', size=3.5) + 
  geom_rect(xmin = 71, xmax=104, ymin=10, ymax=18, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 88, y=15,label= '
Year 1 \n 39 missing outcome \n Year 2 \n 23 missing outcome', size=3.5) + 
  geom_rect(xmin = 71, xmax=104, ymin=-3, ymax=9, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 88, y=4,label= '
Year 1 \n 62 clusters \n 341 children \n Year 2 \n 67 clusters \n 378 children', size=3.5) +
  geom_rect(xmin = 9, xmax=25, ymin=63, ymax=75, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 17, y=70,label= '
Year 1 \n 68 clusters \n 516 children \n Year 2 \n 68 clusters \n 516 children ', size=3.5) +
  annotate('text', x= 3, y=70,label= 'Subsample \n Target', size=3.5) +
  geom_rect(xmin = 6, xmax=28, ymin=32, ymax=62, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 17, y=48,label= '
Year 1 \n 140 children lost to follow-up \n 14 moved \n 16 absent \n 62 withdrew \n 29 no live birth \n 19 child death \n Year 2 \n 0 new children measured  \n 158 children lost to follow-up \n 35 moved \n 3 absent \n 72 withdrew \n 29 no live birth \n 19 child death ', size=3.5) +
  annotate('text', x= 0, y=49,label= 'Follow-up', size=3.3) +
  geom_rect(xmin = 9, xmax=25, ymin=19, ymax=31, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 17, y=26,label= '
Year 1 \n 68 clusters \n 376 children \n Year 2 \n 68 clusters \n 358 children ', size=3.5) +
  annotate('text', x= 2.5, y=26,label= 'Subsample \n Enrollment', size=3.5) +
  geom_rect(xmin = 9, xmax=25, ymin=10, ymax=18, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 17, y=15,label= '
Year 1 \n 57 missing outcome \n Year 2 \n 34 missing outcome', size=3.5) +
  annotate('text', x= 2.5, y=15,label= 'Specimen \n Collection', size=3.5) +
  annotate('text', x= 2.5, y=4,label= 'Analysis', size=3.5) +
  geom_rect(xmin = 9, xmax=25, ymin=-3, ymax=9, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 17, y=4,label= '
Year 1 \n 68 clusters \n 319 children \n Year 2 \n 68 clusters \n 324 children', size=3.5) ->
  p
p

p +
  geom_segment(
    x=50, xend=50, y=100, yend=89.5, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=50, xend=58, y=95, yend=95, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=17, xend=17, y=76, yend=75, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=17, xend=17, y=63, yend=62, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=17, xend=17, y=32, yend=31, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=17, xend=17, y=19, yend=18, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=17, xend=17, y=10, yend=9, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=88, xend=88, y=76, yend=75, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=88, xend=88, y=63, yend=62, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=88, xend=88, y=32, yend=31, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=88, xend=88, y=19, yend=18, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=88, xend=88, y=10, yend=9, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=25, xend=17, y=86.5, yend=86.5, 
    size=0.15, linejoin = "mitre", lineend = "butt") + 
  geom_segment(
    x=17, xend=17, y=86.5, yend=82, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  geom_segment(
    x=75, xend=88, y=86.5, yend=86.5, 
    size=0.15, linejoin = "mitre", lineend = "butt") + 
  geom_segment(
    x=88, xend=88, y=86.5, yend=82, 
    size=0.15, linejoin = "mitre", lineend = "butt",
    arrow = arrow(length = unit(1, "mm"), type= "closed")) ->
  p
p


# YOU MAY NEED TO CHANGE THE FILE PATHS HERE
ggsave(p, file = here("figures/enrollment_figure1.jpg"), height=14, width=9)
