rm(list=ls())

#packages
library(sf)
library(mapplots)
library(scatterpie)
library(maptools)
library(plyr) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggspatial)
library(ggrepel)

# To run this script, change the "mainwd" to wherever this folder
# ("eidolon-dupreanum-henipavirus") is stored on your computer
# Also, make sure to download/clone the "Mada-GIS" folder to 
# your home computer. I recommend putting it in the same parent 
# directory as "eidolon-dupreanum-henipavirus"

# For example, my two folders are stored at:

# "/Users/caraebrook/Documents/R/R_repositories/eidolon-dupreanum-henipavirus/     ...AND
# "/Users/caraebrook/Documents/R/R_repositories/Mada-GIS/

# I keep all my github repos under "R_repositories"

#####################################################################
#####################################################################
# Set wd to data on this computer. Also ID homewd, assuming that 
# Mada-GIS is cloned to the same series of sub-folders
homewd = "/Users/carabrook/Developer/angavokely-virus" 
#should be wherever "eidolon-dupreanum-henipavirus" is stored on your home computer
basewd = paste(strsplit(homewd, "/")[[1]][1:4], collapse = "/")
mapwd = paste0(basewd, "/", "Mada-GIS")
setwd(paste0(homewd, "/", "Fig1/"))



#import madagascar shapfile
name<- paste0(mapwd, "/", "MDG-3/MDG_adm3.shp")
otl_file <- paste(name, sep="") 
orotl_shp <- st_read(otl_file)
#View(orotl_shp)  # Open attribute table
class(orotl_shp)

###import and configuration
# plot mada
# note that this may bog your computer down : I only 
# recommend printing it once to check. If too slow, you can always
# comment out the "print" line and save it temporarily as a pdf instead
# (save script is commented out below the plot)

p1<-ggplot() +  
  geom_sf(color = "lightgoldenrod1", fill = "lightgoldenrod1",data = orotl_shp)+
  coord_sf(xlim = c(42, 60), ylim = c(-26, -11.5), expand = FALSE)+
  theme_bw()+
  theme(plot.margin = unit(c(-1,.5,-1.5,.1),"cm"))+
  xlab("Longitude") + ylab("Latitude") 
#print(p1)
# # 
#   ggsave(file = paste0(homewd, "final-figures/tmp1.pdf"),
#          plot = p1,
#          units="mm",  
#          width=60, 
#          height=55, 
#          scale=3, 
#          dpi=300)
# 



#import henipa data
dat <- read.csv(file = paste0(homewd,"/metadata/NGS_urine_metadat.csv"), header = T, stringsAsFactors = F )
head(dat)

length(unique(dat$sampleid))#206

#only plot urine
dat = subset(dat, sample_type=="urine")# & bat_species=="Eidolon dupreanum") 206 

#and move a few coordinates for visibility
dat$longitude_e[dat$roost_site=="AngavoKely"] <- 47.6
dat$latitude_s[dat$roost_site=="AngavoBe"] <- -19.2


#add age class
#clean class
unique(dat$bat_age_class)

#and rank by rough age
#unique(dat$young_of_year)
dat$age_class <- dat$bat_age_class
dat$age_class[dat$age_class=="P" | dat$age_class=="L"] <- "A"
dat$age_class[dat$age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$age_class[dat$young_of_year=="yes"] <- "J"

#collapse to just one entry per individual (a few bats gave multiple urine specimens)
bat.double = subset(dat, sampleid=="MIZ151" |sampleid=="MIZ275" | sampleid=="MIZ281" |sampleid=="MIZ298" |sampleid=="MIZ300" | sampleid=="MIZ305a" | sampleid=="MIZ305b")
unique(bat.double$henipavirus) #all negative.

doub.list <- dlply(bat.double, .(sampleid))
condense.dat <- function(df){
  df1 <- df[1,]
  return(df1)
  
}
bat.slim <- data.table::rbindlist(lapply(doub.list, condense.dat))


#now collapse those duplicates into one entry for the table
bat.single = subset(dat, sampleid!="MIZ151" &sampleid!="MIZ275" & sampleid!="MIZ281" &sampleid!="MIZ298" &sampleid!="MIZ300" & sampleid!="MIZ305a" & sampleid!="MIZ305b")
dat <- rbind(bat.single, bat.slim)

# now subset the data to just include the columns of interest

map.dat <- dplyr::select(dat,roost_site,latitude_s, longitude_e,
                     collection_date, age_class, bat_sex,
                     bat_species, sampleid, henipavirus)

head(map.dat)
unique(map.dat$roost_site)

#get sites
coordinate <- ddply(map.dat, .(roost_site), summarise, latitude_s=unique(latitude_s), longitude_e=unique(longitude_e))

#now plot the proportion positive at each site

#we will just plot eidolon
#coordinate <-subset(coordinate, roost_site=="AngavoKely" )
coordinate$bat_species <- NA# c("Eidolon dupreanum")
coordinate$bat_species[coordinate$roost_site=="Ambakoana" | coordinate$roost_site=="Makira" | coordinate$roost_site=="Mahialambo" | coordinate$roost_site=="Mahabo" | coordinate$roost_site=="Marovitsika"] <- "Pteropus rufus"
coordinate$bat_species[coordinate$roost_site=="Maromizaha"] <- "Rousettus madagascariensis"
coordinate$bat_species[is.na(coordinate$bat_species)] <- "Eidolon dupreanum"
head(coordinate)

#plot sites on map
# p2<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="#97B5CC",size=1,data=dat)+
#   annotation_scale(location = "bl", width_hint = 0.05) +    # scale
#   annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
#                          pad_x = unit(0.02, "cm"), 
#                          pad_y = unit(0.2, "cm"),        
#                          style = north_arrow_fancy_orienteering)
# #print(p2)
# 
#   ggsave(file = "tmp_map_2.pdf",
#          plot = p2,
#           units="mm",  
#           width=40, 
#           height=60, 
#           scale=3, 
#           dpi=300)
# # 
coordinate$label <- coordinate$bat_species
coordinate$label[coordinate$label=="Rousettus madagascariensis"] <- "Rousettus\nmadagascariensis"
coordinate$label[coordinate$label=="Eidolon dupreanum"] <- "Eidolon\ndupreanum"
coordinate.text = coordinate
coordinate.text$latitude_s = -17.5
coordinate.text$longitude_e = 52
coordinate.text$label = coordinate.text$bat_species





###Grouping data for scatterpie
map.dat$plot_class <- NA
map.dat$plot_class[ dat$henipavirus==1] <- "henipavirus pos"
map.dat$plot_class[ dat$henipavirus==0] <- "henipavirus neg"

pies <- ddply(map.dat, .(bat_species, roost_site, latitude_s, longitude_e, plot_class), summarise, value=length(sampleid))



tot_sum = ddply(pies,.(bat_species, roost_site), summarise,N=sum(value))

pies <- merge(pies, tot_sum, by=c("bat_species", "roost_site"), all.x=T)

pies$plot_class <- factor(pies$plot_class, levels=c( "henipavirus pos", "henipavirus neg"))



###Get the pie data in the right format###
colz = c('henipavirus pos' ="firebrick3", 'henipavirus neg' ="dodgerblue")


# p3<-ggplot() + 
#   geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/500)), 
#                   data = pies, cols="plot_class", long_format=TRUE) +
#   scale_fill_manual(values=colz)
# 
# ggsave(file = "tmp_pies.pdf",
#        plot = p3,
#        units="mm",  
#        width=40, 
#        height=60, 
#        scale=3, 
#        dpi=300)


# 
# # copie of latitude (x.) and longitude (y.)
pies$x2 <- pies$longitude_e
pies$y2 <- pies$latitude_s
# 
# #manually move the pie chart in case there is an overlap (change x and y)


#pies$x2[pies$bat_species== "Eidolon dupreanum"] <- pies$longitude_e[pies$bat_species== "Eidolon dupreanum"] + 4
#pies$y2[pies$bat_species== "Eidolon dupreanum"] <- pies$latitude_s[pies$bat_species== "Eidolon dupreanum"] - 1

pies$x2[pies$roost_site== "AngavoBe"] <- pies$x2[pies$roost_site== "AngavoBe"] + 2.5
pies$y2[pies$roost_site== "AngavoBe"] <- pies$y2[pies$roost_site== "AngavoBe"] - 1.5
pies$x2[pies$roost_site== "AngavoKely"] <- pies$x2[pies$roost_site== "AngavoKely"] + .5
pies$y2[pies$roost_site== "AngavoKely"] <- pies$y2[pies$roost_site== "AngavoKely"] - 3

pies$x2[pies$roost_site== "Ambakoana"] <- pies$x2[pies$roost_site== "Ambakoana"] -2.5
pies$y2[pies$roost_site== "Ambakoana"] <- pies$y2[pies$roost_site== "Ambakoana"] + 2
pies$x2[pies$roost_site== "Mahialambo"] <- pies$x2[pies$roost_site== "Mahialambo"] -2
pies$y2[pies$roost_site== "Mahialambo"] <- pies$y2[pies$roost_site== "Mahialambo"] + 3
pies$x2[pies$roost_site== "Marovitsika"] <- pies$x2[pies$roost_site== "Marovitsika"] -2
pies$y2[pies$roost_site== "Marovitsika"] <- pies$y2[pies$roost_site== "Marovitsika"] + 1
pies$x2[pies$roost_site== "Mahabo"] <- pies$x2[pies$roost_site== "Mahabo"] -1.5
pies$y2[pies$roost_site== "Mahabo"] <- pies$y2[pies$roost_site== "Mahabo"] 
pies$x2[pies$roost_site== "Makira"] <- pies$x2[pies$roost_site== "Makira"] +1.5
pies$y2[pies$roost_site== "Makira"] <- pies$y2[pies$roost_site== "Makira"] -.5
pies$x2[pies$roost_site== "Maromizaha"] <- pies$x2[pies$roost_site== "Maromizaha"] +2
pies$y2[pies$roost_site== "Maromizaha"] <- pies$y2[pies$roost_site== "Maromizaha"] +1


head(pies)



#plot pie chart 
#loko<-c("Rousettus madagascariensis"="#B200ED","Eidolon dupreanum"="#7FFF00","Pteropus rufus"="#0000FF")

#this is Fig1A

shapez = c('Pteropus rufus' = 21, 'Eidolon dupreanum' = 24, 'Rousettus madagascariensis' = 22)
colz1 = c("Pteropus rufus" = "magenta", "Eidolon dupreanum" = "mediumseagreen", "Rousettus madagascariensis" = "darkorange1")

colz = c('henipavirus pos' ="firebrick3", 'henipavirus neg' ="dodgerblue")


Fig1a <- p1+
  annotate("segment", x=pies$longitude_e, xend=pies$x2,y=pies$latitude_s,yend=pies$y2,size=.5, linetype=5, color="black")+ # put the lines
  geom_point(aes(x=longitude_e, y=latitude_s, fill=bat_species, shape = bat_species), color="black", stroke=1,
             size=3.5,data=coordinate)+
  scale_shape_manual(values=shapez, guide=guide_legend(label.theme=element_text(face="italic", angle = 0, size=10))) +
  scale_fill_manual(values=colz1, guide=guide_legend(label.theme=element_text(face="italic", angle = 0, size=10))) +
  annotation_scale(location = "bl", width_hint = 0.05) +    #scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.03, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)+
  geom_text_repel(segment.colour="black")+
  ggnewscale::new_scale_fill() +
  geom_scatterpie(aes(x=x2, y=y2, r=(N/50)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(0.7,.5,.3,.1),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.8,.87),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size=10)) +
  scale_fill_manual(values=colz) +
  coord_sf(xlim=c(43,56))+
  geom_scatterpie_legend((c(50,100)/50),
                         x=52, y=-24, 
                         n=2,
                         labeller = function(x) paste((x)*50,"indiv"))



ggsave(file = paste0(homewd, "/final-figures/Fig1.pdf"),
       plot=Fig1a,
       units="mm",  
       width=55, 
       height=65, 
       scale=3, 
       dpi=300)


#other file types

ggsave(file = paste0(homewd, "/final-figures/Fig1.tiff"),
       plot=Fig1a,
       units="mm",  
       width=55, 
       height=65, 
       scale=3, 
       dpi=300)
