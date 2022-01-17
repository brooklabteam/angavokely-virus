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
setwd(paste0(homewd, "/Fig2"))

#and load the sequences
seq <- read.fasta(file ="Orthoparamyxovirinae.fasta", forceDNAtolower = F, as.string = T)
dat <- data.frame(names=names(seq))
dat$n <- sapply(strsplit(dat$names, "_"),length)

#and convert to just the accession number
dat$new_name <- NA
dat$new_name[dat$n==1] <- sapply(strsplit(dat$names[dat$n==1], "_"), function(x) x[[1]])
dat$new_name[dat$n>1] <- paste0(sapply(strsplit(dat$names[dat$n>1], "_"), function(x) x[[1]]), "_", sapply(strsplit(dat$names[dat$n>1], "_"), function(x) x[[2]]))
dat$new_name[dat$new_name=="408_forced"] <- "AngV"
dat$new_name <-sub(pattern=".", replacement = "_", x=dat$new_name, fixed=T)


write.fasta(seq, names = dat$new_name, file.out = "Orthoparamyxovirinae_newnames.fasta")

#now send to alignment, then Modeltest-NG and RAxML