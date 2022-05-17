rm(list=ls())


library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(lubridate)
#library(rBt)
library(treeio)


#make Bayesian timetree from Henipavirus strict molecular clock model

#first, read in the tree
homewd= "/Users/carabrook/Developer/angavokely-virus"
setwd(paste0(homewd, "/Fig4"))

tree <- read.beast(file = paste0(homewd, "/Fig4/henipaAVG.tree"))

treedat <- cbind.data.frame(tip_name = tree@phylo$tip.label)
treedat$beast_name <-treedat$tip_name

#tree$node.label <- round(tree$posterior,2)
#treedat <- cbind.data.frame(tip_name = tree$tip.label)
treedat$Accession_Num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])
treedat$Accession_Num[treedat$Accession_Num=="NC"] <- c("NC_009021", "NC_048212")
#names(treedat)[names(treedat)=="tip_name"] <- "beast_name"


#and load data of corresponding tree
dat <- read.csv(file = paste0(homewd, "/Fig4/beast_seq.csv"), header = T, stringsAsFactors = F)
head(dat)

#and merge on tiplabel
tree@phylo$tip.label[tree@phylo$tip.label=="408_forced_SPADES_contig_16740"] <- "AngV"
treedat$tip_name[treedat$tip_name=="408_forced_SPADES_contig_16740"] <- "AngV"
treedat$beast_name[treedat$beast_name=="408_forced_SPADES_contig_16740"] <- "AngV"
treedat$Accession_Num[treedat$beast_name=="AngV"] <- "AngV"

mergedat <- merge(treedat, dat, by="Accession_Num", all.x = T, sort = F)
head(mergedat)
#mergedat$Collection_Date <- as.Date(mergedat$Collection_Date, format = "%d-%b-%Y")
mergedat$Collection_Date <- as.Date(mergedat$Collection_Date)
mergedat$year <- year(mergedat$Collection_Date)
#and add a new label
mergedat$new_label <- paste0(mergedat$Accession_Num, " | ", mergedat$Name, " | ", mergedat$Locality,
                             " | ",  mergedat$year, " | ", mergedat$Host)


mrsd.dat <- as.Date(max(mergedat$Collection_Date))

#convert x to date in the same form - "0" should equal the most recent date
p = ggtree(tree) + theme_tree2()
p$data$x = p$data$x - max(p$data$x) #+ 463

p1 <- p + theme_tree2() + 
    scale_x_continuous(breaks=c(-10000, -8000, -6000, -4000, -2000, 0),
    labels=c(10000, 8000, 6000, 4000, 2000,  0)) +
    xlab("years to MRCA")

names(mergedat)[names(mergedat)=="Name"] <- "virus"


tree@data$posterior[tree@data$posterior<0.9] <- 0
tree@data$posterior[tree@data$posterior>=0.9] <- 1
tree@data$posterior <- round(tree@data$posterior, 2)
tree@data$posterior <- as.character(tree@data$posterior)
tree@data$posterior[tree@data$posterior=="0"] <- "<0.9"
tree@data$posterior[tree@data$posterior=="1"] <- ">0.9"
tree@data$posterior <- factor(tree@data$posterior, levels= c("<0.9", ">0.9"))

tree@phylo$tip.label <- mergedat$new_label
mergedat$new_label
mergedat <- dplyr::select(mergedat, new_label, names(mergedat)[1:9])

p = ggtree(tree) + theme_tree2()
p$data$x = p$data$x - max(p$data$x) #+ 463
#remove HPD that are the more recent -- only include from >500 yrs
#p$data$height_0.95_HPD[p$data$height<400] <- NA
#add the big one
#p$data$length_0.95_HPD[p$data$height==max(p$data$height)] <- p$data$height_range[p$data$height==max(p$data$height)] 

p1 <- p + theme_tree2() + 
  scale_x_continuous(breaks=c(-10000, -8000, -6000, -4000, -2000, 0),
                     labels=c(10000, 8000, 6000, 4000, 2000,  0)) +
  xlab("years to MRCA")


posfilz <- c('<0.9'="white", '>0.9'="black")

mergedat$novel = "no"
mergedat$novel[mergedat$Locality=="Madagascar"] <- "yes"

colz2 = c('yes' =  "yellow", 'no' = "white")



p2 <- p1 %<+% mergedat + 
  #geom_tiplab(size=2, nudge_x = 500) + #geom_nodelab(size=2,nudge_x = -15, nudge_y = .7) +
  theme(legend.position = c(.15,.75)) +
  geom_range(range='height_0.95_HPD', color='red', alpha=.5, size=1) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  geom_tippoint(aes(color=virus)) +
  scale_fill_manual(values=posfilz) +
  guides( fill_continuous = guide_legend(order = 2, nrow = 1),col = guide_legend(order = 1)) +
  coord_cartesian(clip = "off", xlim=c(-15000, 9000)) +ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=2, nudge_x=500) + scale_fill_manual(values=colz2)


  


ggsave(file = paste0(homewd, "/final-figures/SuppFig1.png"),
       plot=p2,
       units="mm",  
       width=70, 
       height=90, 
       #limitsize = F,
       scale=3)#, 

#and save this as the Supp plot


#and collapse some clades

#nodeNiV <- MRCA(tree, which(tree@phylo$tip.label == "MK673565 | NiV | Bangladesh | 2004 | Homo_sapiens"),which(tree@phylo$tip.label == "MK801755 | NiV | Cambodia | 2003 | Pteropus_lylei"))
nodeNiVbat <- MRCA(tree, which(tree@phylo$tip.label == "AJ564623 | NiV | Malaysia | 1999 | Sus_scrofa_domesticus"),which(tree@phylo$tip.label == "AJ627196 | NiV | Malaysia | 1999 | Sus_scrofa_domesticus"))
nodeNiVIndia <- MRCA(tree, which(tree@phylo$tip.label == "MH523642 | NiV | India | 2018 | Homo_sapiens"),which(tree@phylo$tip.label == "MN549409 | NiV | India | 2019 | Pteropus_medius"))
nodeNiVBangBat <- MRCA(tree, which(tree@phylo$tip.label == "MK575066 | NiV | Bangladesh | 2013 | Pteropus_medius"),which(tree@phylo$tip.label == "MK575060 | NiV | Bangladesh | 2013 | Pteropus_medius"))
nodeNiVBangHum <- MRCA(tree, which(tree@phylo$tip.label == "MK673590 | NiV | Bangladesh | 2014 | Homo_sapiens"),which(tree@phylo$tip.label == "MK673565 | NiV | Bangladesh | 2004 | Homo_sapiens"))
nodeNiVBangHum2 <- MRCA(tree, which(tree@phylo$tip.label == "MK673586 | NiV | Bangladesh | 2015 | Homo_sapiens"),which(tree@phylo$tip.label == "MK673579 | NiV | Bangladesh | 2012 | Homo_sapiens"))

#nodeHeVhorse <- MRCA(tree, which(tree@phylo$tip.label == "HM044319 | HeV | Australia | 2007 | Equus_f_caballus"),which(tree@phylo$tip.label == "HM044320 | HeV | Australia | 2008 | Equus_f_caballus"))
nodeHeV2007 <- MRCA(tree, which(tree@phylo$tip.label == "JN255804 | HeV | Australia | 2007 | Pteropus_alecto/Pteropus_poliocephalus"),which(tree@phylo$tip.label == "HM044321 | HeV | Australia | 2007 | Equus_f_caballus"))
nodeHeV1994 <- MRCA(tree, which(tree@phylo$tip.label == "MN062017 | HeV | Australia | 1994 | Equus_f_caballus"),which(tree@phylo$tip.label == "AF017149 | HeV | Australia | 1994 | Homo_sapiens"))
nodeHeV2009 <- MRCA(tree, which(tree@phylo$tip.label == "JN255803 | HeV | Australia | 2009 | Pteropus_alecto/Pteropus_poliocephalus"),which(tree@phylo$tip.label == "JN255800 | HeV | Australia | 2009 | Pteropus_alecto"))


#highlight Mada

#make tip labels a little bigger

p2 <- p1 %<+% mergedat + 
  #geom_tiplab(size=3, nudge_x = 500) + #geom_nodelab(size=2,nudge_x = -15, nudge_y = .7) +
  theme(legend.position = c(.1,.75)) +
  geom_range(range='height_0.95_HPD', color='red', alpha=.5, size=1) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=2, stroke=.1) +
  geom_tippoint(aes(color=virus), size=3) +
  scale_fill_manual(values=posfilz) +
  guides( fill_continuous = guide_legend(order = 2, nrow = 1),col = guide_legend(order = 1)) +
  coord_cartesian(clip = "off", xlim=c(-14500, 20000)) +ggnewscale::new_scale_fill() +
  geom_tiplab(aes(fill = novel), geom = "label", label.size = 0, 
              alpha=.3,  show.legend=F, size=3.5, nudge_x=500) + scale_fill_manual(values=colz2)





p2.1 <- collapse(p2, node=nodeNiVbat) #this node is malaysia 1999 
p2.2 <- collapse(p2.1, node=nodeNiVIndia) #this is India 2018/2019
p2.3 <- collapse(p2.2, node=nodeNiVBangBat) #this is Bangladesh bats 2013
p2.4 <- collapse(p2.3, node=nodeNiVBangHum) #this is Bangladesh humans 
p2.5 <- collapse(p2.4, node=nodeNiVBangHum2) #this is Bangladesh humans 2
p2.6 <- collapse(p2.5, node=nodeHeV2007) #this is 2007 horse and bat
p2.7 <- collapse(p2.6, node=nodeHeV1994) # this is 1997 horse and human
p2.8 <- collapse(p2.7, node=nodeHeV2009) # this is 2009 just bat

#this is the NiV color
scales::hue_pal()(8)[8]

#and the HeV color
scales::hue_pal()(8)[6]


#now add node labels for the clusters and color these nodes pink with s new symbol
p2.9 <- p2.8 + 
        geom_point2(aes(subset=(node==nodeNiVbat)), shape=17, size=3, fill="#FF61CC", color="#FF61CC")+
  geom_point2(aes(subset=(node==nodeNiVIndia)), shape=17, size=3, fill="#FF61CC", color="#FF61CC") +
  geom_point2(aes(subset=(node==nodeNiVBangBat)), shape=17, size=3, fill="#FF61CC", color="#FF61CC") +
  geom_point2(aes(subset=(node==nodeNiVBangHum)), shape=17, size=3, fill="#FF61CC", color="#FF61CC") +
  geom_point2(aes(subset=(node==nodeNiVBangHum2)), shape=17, size=3, fill="#FF61CC", color="#FF61CC") +
  geom_point2(aes(subset=(node==nodeHeV2007)), shape=17, size=3, fill="#00A9FF", color="#00A9FF") +
  geom_point2(aes(subset=(node==nodeHeV1994)), shape=17, size=3, fill="#00A9FF", color="#00A9FF") +
  geom_point2(aes(subset=(node==nodeHeV2009)), shape=17, size=3, fill="#00A9FF", color="#00A9FF")



p3 <- p2.9 + geom_nodelab(aes(subset=(node==nodeNiVbat), label = "NiV | Malaysia | 1999 | Sus_scrofa/Canis_familiaris/Homo_sapiens"), geom="text", size=3.5, nudge_x=10500)+
               geom_nodelab(aes(subset=(node==nodeNiVIndia), label = "NiV | India | 2018-2019 | Homo_sapiens/Pteropus_medius"), geom="text", size=3.5, nudge_x=9200)+
               geom_nodelab(aes(subset=(node==nodeNiVBangBat), label = "NiV | Bangladesh | 2013 | Pteropus_medius"), geom="text", size=3.5, nudge_x=8500)+
               geom_nodelab(aes(subset=(node==nodeNiVBangHum), label = "NiV | Bangladesh | 2004-2014 | Homo_sapiens"), geom="text", size=3.5, nudge_x=7500)+
               geom_nodelab(aes(subset=(node==nodeNiVBangHum2), label = "NiV | Bangladesh/India | 2004-2015 | Homo_sapiens"), geom="text", size=3.5, nudge_x=8500) +
               geom_nodelab(aes(subset=(node==nodeHeV2007), label = "HeV | Australia | 2007 | Equus_f_caballus/Pteropus_sp"), geom="text", size=3.5, nudge_x=8500)+
               geom_nodelab(aes(subset=(node==nodeHeV1994), label = "HeV | Australia | 1994 | Equus_f_caballus/Homo_sapiens"), geom="text", size=3.5, nudge_x=9000)+
               geom_nodelab(aes(subset=(node==nodeHeV2009), label = "HeV | Australia | 2009 | Pteropus_sp"), geom="text", size=3.5, nudge_x=6000)


ggsave(file = paste0(homewd, "/final-figures/Fig4B.png"),
       plot=p3,
       units="mm",  
       width=70, 
       height=80, 
       #limitsize = F,
       scale=3)#, 

#and to find the MRCA
#divergence of AngV:

nodebranchAngV = MRCA(tree, which(tree@phylo$tip.label == "AngV | AngV | Madagascar | 2019 | Eidolon_dupreanum"),which(tree@phylo$tip.label == "HQ660129 | GhV | Ghana | 2009 | Eidolon_helvum"))
nodeBasal = MRCA(tree, which(tree@phylo$tip.label == "AngV | AngV | Madagascar | 2019 | Eidolon_dupreanum"),which(tree@phylo$tip.label == "KF278639 | MojV | China | 2012 | Rattus_rattus"))

dat.phylo <- p3$data
head(dat.phylo)

#and get the node of interest:
subset(dat.phylo, node==nodebranchAngV) #9794 years
subset(dat.phylo, node==nodeBasal) #11195 years

#and the 95% HPD
dat.phylo$height_0.95_HPD[[nodebranchAngV]] #6519.024 14023.649
dat.phylo$height_0.95_HPD[[nodeBasal]] #7350.913 15904.564
