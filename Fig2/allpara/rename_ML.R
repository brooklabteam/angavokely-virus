rm(list=ls())


library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
library(seqinr)


#make Bayesian timetree from Henipavirus strict molecular clock model

#first, read in the tree
homewd= "/Users/carabrook/Developer/angavokely-virus"
setwd(paste0(homewd, "/Fig2/allpara/"))

#and load the sequences
seq <- read.fasta(file ="allpara.fasta", forceDNAtolower = F, as.string = T)

dat <- data.frame(names=c(unlist(lapply(seq, function(x) attr(x,"Annot")))))
dat$names <-sub(pattern=">", replacement = "", x=dat$names, fixed=T)
dat$names <-sub(pattern=".", replacement = "_", x=dat$names, fixed=T)
dat$n <- sapply(strsplit(dat$names, "_"),length)
#and convert to just the accession number
dat$new_name <- NA
dat$new_name[dat$n==1] <- sapply(strsplit(dat$names[dat$n==1], "_"), function(x) x[[1]])

dat$new_name[dat$n>1] <- paste0(sapply(strsplit(dat$names[dat$n>1], "_"), function(x) x[[1]]), "_", sapply(strsplit(dat$names[dat$n>1], "_"), function(x) x[[2]]))
dat$word_name <- NA
dat$word_name[dat$n>2] <- sapply(strsplit(dat$names[dat$n>2], "_"), function(x) x[[3]])
dat$word_name <- sub(pattern="1 |", replacement = "", x=dat$word_name, fixed=T)
dat$word_name <- sub(pattern="2 |", replacement = "", x=dat$word_name, fixed=T)
dat$word_name <- sub(pattern="3 |", replacement = "", x=dat$word_name, fixed=T)
write.csv(dat, file = "ML_dat_allpara.csv", row.names = F)

#and load
dat <- read.csv(file =  "ML_dat_allpara.csv", stringsAsFactors = F, header = T)
head(dat)
dat$add_name <- paste0(dat$new_name, "_", dat$word_name)
dat$add_name[dat$add_name=="AngV_AngV"] <- "AngV"

write.fasta(seq, names = dat$add_name, file.out = "AllParamyxovirus_newnames.fasta")

#now send to alignment, then Modeltest-NG and RAxML