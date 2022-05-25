rm(list = ls())

library(plyr)
library(dplyr)

homewd = "/Users/carabrook/Developer/angavokely-virus"
setwd(homewd)
ang.dat <- read.csv(paste0(homewd, "/metadata/czb_meta_bat_urine.csv"))
head(ang.dat)
length(unique(ang.dat$czb_id)) #213
length(unique(ang.dat$sampleID)) #205

dat.list <- dlply(ang.dat, .(sampleID))
get.samples <- function(df){
  df.out <- cbind.data.frame(sampleid = unique(df$sampleID), num_samples = nrow(df))
  return(df.out)
}

dat.new <- data.table::rbindlist(lapply(dat.list, get.samples))
dat.new
subset(dat.new, num_samples>1)
