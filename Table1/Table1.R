rm(list = ls())

library(plyr)
library(dplyr)

homewd = "/Users/carabrook/Developer/angavokely-virus"
setwd(homewd)
ang.dat <- read.csv(paste0(homewd, "/metadata/NGS_urine_metadat.csv"))
head(ang.dat)
length(unique(ang.dat$sampleid)) #206 unique bats
length(unique(ang.dat$czb_id)) #213


dat.list <- dlply(ang.dat, .(sampleid))
get.samples <- function(df){
  df.out <- cbind.data.frame(sampleid = unique(df$sampleid), num_samples = nrow(df))
  return(df.out)
}

dat.new <- data.table::rbindlist(lapply(dat.list, get.samples))
dat.new
subset(dat.new, num_samples>1) #7 bats sampled twice

#for those bats sampled twice, do they ever have different outcomes for the different samples?
#should be no because they are all rousettus, but let's just double check
bat.double = subset(ang.dat, sampleid=="MIZ151" |sampleid=="MIZ275" | sampleid=="MIZ281" |sampleid=="MIZ298" |sampleid=="MIZ300" | sampleid=="MIZ305a" | sampleid=="MIZ305b")
unique(bat.double$henipavirus) #all negative.

doub.list <- dlply(bat.double, .(sampleid))
condense.dat <- function(df){
  df1 <- df[1,]
  return(df1)
  
}
bat.slim <- data.table::rbindlist(lapply(doub.list, condense.dat))


#now collapse those duplicates into one entry for the table
bat.single = subset(ang.dat, sampleid!="MIZ151" &sampleid!="MIZ275" & sampleid!="MIZ281" &sampleid!="MIZ298" &sampleid!="MIZ300" & sampleid!="MIZ305a" & sampleid!="MIZ305b")
bat.table <- rbind(bat.single, bat.slim)
head(bat.table)
bat.table$collection_date <- as.Date(bat.table$collection_date, format = "%m/%d/%y")
library(lubridate)
bat.table$year <- year(bat.table$collection_date)
bat.table$month <- month(bat.table$collection_date)
unique(bat.table$roost_site)
unique(bat.table$bat_species)
unique(bat.table$bat_sex)
#and split by site, species, season, sex
table1.start <- ddply(bat.table, .(bat_species, roost_site, year, month, bat_sex), summarise, N_sampled = length(henipavirus), N_positive=sum(henipavirus))
#and combine a few into the same season (e.g. wet season spans beginning of one year to start of next)
table1.start$season <- NA
table1.start$season[table1.start$month>=6 & table1.start$month <= 9] <- "dry"
table1.start$season[table1.start$month > 9] <- "wet-late"
table1.start$season[table1.start$month < 6] <- "wet-early"
unique(table1.start$season)
table1.start$season_sum <- NA 
table1.start$season_sum[table1.start$season=="dry"] <- paste0(table1.start$year[table1.start$season=="dry"], "-", table1.start$season[table1.start$season=="dry"] )
subset(table1.start, is.na(season_sum))
table1.start$season_sum[table1.start$season=="wet-late" & table1.start$year==2013] <- "2013-wet"
table1.start$season_sum[table1.start$season=="wet-late" & table1.start$year==2014 | table1.start$season=="wet-early" & table1.start$year==2015] <- "'14/'15-wet"
table1.start$season_sum[table1.start$season=="wet-late" & table1.start$year==2015 | table1.start$season=="wet-early" & table1.start$year==2016] <- "'15/'16-wet"
table1.start$season_sum[table1.start$season=="wet-late" & table1.start$year==2016 | table1.start$season=="wet-early" & table1.start$year==2017] <- "'17/'17-wet"
table1.start$season_sum[table1.start$season=="wet-late" & table1.start$year==2017 | table1.start$season=="wet-early" & table1.start$year==2018] <- "'17/'18-wet"
table1.start$season_sum[table1.start$season=="wet-late" & table1.start$year==2018 | table1.start$season=="wet-early" & table1.start$year==2019] <- "'18/'19-wet"

#can any be combined?
unique(table1.start$season_sum)

#and combine
table1 <- ddply(table1.start, .(bat_species, roost_site, season_sum, bat_sex), summarise, N_sampled = sum(N_sampled), N_positive=sum(N_positive))

#and write to files
write.csv(table1, file = paste0(homewd, "/final-figures/table1.csv"), row.names = F)

#and open in excel for some manual edits too
