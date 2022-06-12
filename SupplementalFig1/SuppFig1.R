rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(seqinr)
library(ggmsa)


#load the msa of just the henipavirus W genoms

homewd = "/Users/carabrook/Developer/angavokely-virus/"
setwd(paste0(homewd, "SupplementalFig1")) 

#here's the whole protein
Wmsa = tidy_msa(msa="WproteinAlignment.fasta")
head(Wmsa)

#and here's just the C-terminal
CtermWmsa = tidy_msa(msa="WproteinCterminalAlignment.fasta")
head(CtermWmsa)

#assign new names to each
Wmsa$name <- as.character(Wmsa$name)
unique(Wmsa$name)
Wmsa$name[Wmsa$name=="AngV_addGG_W_translation"] <- "AngV"
Wmsa$name[Wmsa$name=="NiV_MW535746_W_Translation Nipah henipavirus isolate NiV/TH/P.lylei/2017/B17640.GUL, complete genome"] <- "NiV"
Wmsa$name[Wmsa$name=="HeV_MZ229748_W_translation Hendra henipavirus isolate HeV-G2/Australia/Pteropus poliocephalus/2020-03, complete genome"] <- "HeV"
Wmsa$name[Wmsa$name=="MojV_KF278639_addGG_W_translation Mojiang virus isolate Tongguan1, complete genome"] <- "MojV"
Wmsa$name[Wmsa$name=="GhV_HQ660129_addGG_W_translation"] <- "GhV"
Wmsa$name[Wmsa$name=="GAKV_MZ574407-addGG-W_Translation UNVERIFIED: Henipavirus sp. strain Cl17-32, complete genome"] <- "GAKV"
Wmsa$name[Wmsa$name=="DARV_MZ574409-addGG-W_Translation Henipavirus sp. strain Cs17-65, complete genome"] <- "DARV"
unique(Wmsa$name)
names(Wmsa)[names(Wmsa)=="character"] <- "aa"



#assign new names to each
CtermWmsa$name <- as.character(CtermWmsa$name)
unique(CtermWmsa$name)
CtermWmsa$name[CtermWmsa$name=="AngV_addGG_W_translation_extraction"] <- "AngV"
CtermWmsa$name[CtermWmsa$name=="NiV_MW535746_W_Translation_extraction Nipah henipavirus isolate NiV/TH/P.lylei/2017/B17640.GUL, complete genome"] <- "NiV"
CtermWmsa$name[CtermWmsa$name=="HeV_MZ229748_W_translation_extraction Hendra henipavirus isolate HeV-G2/Australia/Pteropus poliocephalus/2020-03, complete genome"] <- "HeV"
CtermWmsa$name[CtermWmsa$name=="MojV_KF278639_addGG_W_translation_extraction Mojiang virus isolate Tongguan1, complete genome"] <- "MojV"
CtermWmsa$name[CtermWmsa$name=="GhV_HQ660129_addGG_W_translation_extraction"] <- "GhV"
CtermWmsa$name[CtermWmsa$name=="GAKV_MZ574407-addGG-W_Translation_extraction UNVERIFIED: Henipavirus sp. strain Cl17-32, complete genome"] <- "GAKV"
CtermWmsa$name[CtermWmsa$name=="DARV_MZ574409-addGG-W_Translation_extraction Henipavirus sp. strain Cs17-65, complete genome"] <- "DARV"
unique(CtermWmsa$name)
names(CtermWmsa)[names(CtermWmsa)=="character"] <- "aa"


head(Wmsa)
unique(Wmsa$aa)


#msa.dat$clade[msa.dat$clade=="African Eidolon"] <- 
#colz = c("A" = "seagreen3", "T" = "cornflowerblue", "G"="tomato", "C" = "plum3", "-" = "white")# "N"= "mediumpurple1", "p10"="red", "NS7" = "magenta")

whole.msa.plot <- ggplot(data=Wmsa) + 
  geom_tile(aes(x=position, y=name, fill=aa, color=aa),width=0.9, height=0.9) +
  geom_text(aes(x=position, y=name, label=aa))+
  theme_bw() + 
  #scale_fill_manual(values=colz) +
  #scale_color_manual(values=colz) +
  theme(axis.title.y = element_blank(), 
        #panel.background = element_rect(fill="black"),
        #panel.border = element_line(color="black", size=1),
        panel.grid = element_blank(),
        axis.text = element_text(size=14),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=16),
        legend.text = element_text(size=10),
        plot.margin = unit(c(.1, 2,.1, .1), "cm"),
        legend.position = "bottom", legend.direction = "horizontal",legend.title=element_blank()) +
  xlab("aa position")

colz = scales::hue_pal()(length(unique(CtermWmsa$aa)))
names(colz) <- unique(CtermWmsa$aa)
names(colz)[names(colz)=="-"] <- "white"

#and control the order to match that from Fig 3
CtermWmsa$name <- factor(CtermWmsa$name, levels = rev(c("AngV", "NiV", "HeV", "MojV", "GhV","GAKV", "DARV")))
msa.plot <- ggplot(data=CtermWmsa) + 
  geom_tile(aes(x=position, y=name, fill=aa, color=aa),width=0.9, height=0.9, show.legend = F) +
  geom_text(aes(x=position, y=name, label=aa), show.legend = F)+
  theme_bw() + 
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz) +
  theme(axis.title.y = element_blank(), 
        #panel.background = element_rect(fill="black"),
        #panel.border = element_line(color="black", size=1),
        panel.grid = element_blank(),
        axis.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        legend.text = element_text(size=10),
        plot.margin = unit(c(.1, 2,.1, .1), "cm"),
        legend.position = "bottom", legend.direction = "horizontal",legend.title=element_blank()) +
  xlab("amino acid position\n(relative to AngV)") + scale_x_continuous(breaks = c(1,10,26,47,68,71), labels = c(378,387,403,404, 409, 412), position = "top")


ggsave(file = paste0(homewd, "/final-figures/SuppFig1.png"),
       plot = msa.plot,
       units="mm",  
       width=100, 
       height=60, 
       #limitsize = F,
       scale=3)#, 

