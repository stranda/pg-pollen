bv=readRDS("output/bvs_n200_v4.1.RDS")

library(dplyr)
library(ggplot2)

q2.5 <- function(x) quantile(x,0.025)
q97.5 <- function(x) quantile(x,0.975)
bvsum <- bv %>% select(timeFrom, timeTo, timeSpan, centroidVelocity, nsQuantVelocity_quant0p05,nsQuantVelocity_quant0p95) %>%group_by(timeFrom, timeTo, timeSpan) %>% summarise_all(list(mean,q2.5,q97.5))

names(bvsum) = gsub("fn1","mn",names(bvsum))
names(bvsum) = gsub("fn2","q2.5",names(bvsum))
names(bvsum) = gsub("fn3","q97.5",names(bvsum))
write.table(file="pollenBVs.csv",sep=",",row.names=F,bvsum)

ggplot(bvsum,aes(x=timeFrom, y=nsQuantVelocity_quant0p05_mn)) +geom_point() 
