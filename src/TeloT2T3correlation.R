#Correlation of TS2 and TS3

d<- readRDS("bangladesh-cleaned-master-data.rds")

t2t3cor <- cor.test(d$TS_t3_Z, d$TS_t2_Z)


t2t3cor
