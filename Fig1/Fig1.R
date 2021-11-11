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
homewd = "/Users/caraebrook/Documents/R/R_repositories/eidolon-dupreanum-henipavirus/" 
#should be wherever "eidolon-dupreanum-henipavirus" is stored on your home computer
basewd = paste(strsplit(homewd, "/")[[1]][1:6], collapse = "/")
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



#import CoV data
dat <- read.csv(file = paste0(homewd,"/metadata/NGS_urine_metadat.csv"), header = T, stringsAsFactors = F )
head(dat)

#only plot urine
dat = subset(dat, sample_type=="urine" & bat_species=="Eidolon dupreanum")

#add age class
#clean class
unique(dat$bat_age_class)

#and rank by rough age
unique(dat$young_of_year)
dat$age_class <- dat$bat_age_class
dat$age_class[dat$age_class=="P" | dat$age_class=="L"] <- "A"
dat$age_class[dat$age_class=="NL" | dat$young_of_year=="no"] <- "A"
dat$age_class[dat$young_of_year=="yes"] <- "J"

# now subset the data to just include the columns of interest

map.dat <- dplyr::select(dat,roost_site,latitude_s, longitude_e,
                     collection_date, age_class, bat_sex,
                     bat_species, sampleid, henipavirus)

head(map.dat)
unique(map.dat$roost_site)

#get sites
coordinate <- ddply(map.dat, .(roost_site), summarise, latitude_s=unique(latitude_s), longitude_e=unique(longitude_e))
#we will just plot eidolon
coordinate <-subset(coordinate, roost_site=="AngavoKely" )
coordinate$bat_species <- c("Eidolon dupreanum")
head(coordinate)

#plot sites on map
p2<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="#97B5CC",size=1,data=dat)+
  annotation_scale(location = "bl", width_hint = 0.05) +    # scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.02, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)
#print(p2)

#  ggsave(file = "tmp_map_2.pdf",
#         plot = p2,
#          units="mm",  
#          width=40, 
#          height=60, 
#          scale=3, 
#          dpi=300)
# # 
coordinate$label <- coordinate$bat_species
#coordinate$label[coordinate$label=="Rousettus madagascariensis"] <- "Rousettus\nmadagascariensis"
coordinate$label[coordinate$label=="Eidolon dupreanum"] <- "Eidolon\ndupreanum"
coordinate.text = coordinate
coordinate.text$latitude_s = -17.5
coordinate.text$longitude_e = 52
coordinate.text$label = coordinate.text$bat_species
#load GPS point and label
p2b<-p1+geom_point(aes(x=longitude_e, y=latitude_s),color="black",size=1,data=coordinate)+
  geom_text(data= coordinate.text,                       #### Labeling
            aes(x=longitude_e, y=latitude_s, label=label),
            fontface="italic",
            color = "#1B262C", size=4,
            #nudge_x = c(-2,-3.6,8),
            #nudge_y = c(3,-1.1,-.3),
            check_overlap = T)+
  annotation_scale(location = "bl", width_hint = 0.05) +    #scale
  annotation_north_arrow(location = "tl", which_north = "true",#north arrow     
                         pad_x = unit(0.03, "cm"), 
                         pad_y = unit(0.2, "cm"),        
                         style = north_arrow_fancy_orienteering)+
  geom_text_repel(segment.colour="black")+
  theme_bw() +theme(panel.grid = element_blank(), 
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(-1,.5,-1.5,.1),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.26,.90),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.background = element_rect(color="gray",size = .1),
                    legend.text = element_text(size = 9,face = "italic"))
#print(p2b)
# # # 
#     ggsave(file = "tmp_map_2b.pdf",
#            plot = p2b,
#            units="mm",  
#            width=40, 
#            height=60, 
#            scale=3, 
#            dpi=300)
# # # 


#plot one site per species
#dat$roost_site[dat$bat_species=="Pteropus rufus" & dat$roost_site=="Marovitsika" |dat$bat_species=="Pteropus rufus" & dat$roost_site=="Mahialambo"] <- "Ambakoana"
map.dat$roost_site[map.dat$bat_species=="Eidolon dupreanum"] <- "AngavoKely"
#dat$roost_site[dat$bat_species=="Rousettus madagascariensis"] <- "Maromizaha"

#dat$longitude_e[dat$roost_site=="Ambakoana"] <- coordinate$longitude_e[coordinate$roost_site=="Ambakoana"]
#dat$longitude_e[dat$roost_site=="Mahabo"] <- coordinate$longitude_e[coordinate$roost_site=="Mahabo"]
#dat$longitude_e[dat$roost_site=="Makira"] <- coordinate$longitude_e[coordinate$roost_site=="Makira"]
map.dat$longitude_e[map.dat$roost_site=="AngavoKely"] <- coordinate$longitude_e[coordinate$roost_site=="AngavoKely"]
#dat$longitude_e[dat$roost_site=="Maromizaha"] <- coordinate$longitude_e[coordinate$roost_site=="Maromizaha"]


#dat$latitude_s[dat$roost_site=="Ambakoana"] <- coordinate$latitude_s[coordinate$roost_site=="Ambakoana"]
map.dat$latitude_s[map.dat$roost_site=="AngavoKely"] <- coordinate$latitude_s[coordinate$roost_site=="AngavoKely"]
#dat$latitude_s[dat$roost_site=="Maromizaha"] <- coordinate$latitude_s[coordinate$roost_site=="Maromizaha"]


###Grouping data for scatterpie
map.dat$plot_class <- NA
map.dat$plot_class[ dat$henipavirus==1] <- "henipavirus pos"
map.dat$plot_class[ dat$henipavirus==0] <- "henipavirus neg"

pies <- ddply(map.dat, .(bat_species, roost_site, latitude_s, longitude_e, plot_class), summarise, value=length(sampleid))



tot_sum = ddply(pies,.(bat_species), summarise,N=sum(value))

pies <- merge(pies, tot_sum, by=c("bat_species"), all.x=T)

pies$plot_class <- factor(pies$plot_class, levels=c( "henipavirus pos", "henipavirus neg"))



###Get the pie data in the right format###
colz = c('henipavirus pos' ="firebrick4", 'henipavirus neg' ="dodgerblue")


p3<-ggplot() + 
  geom_scatterpie(aes(x=longitude_e, y=latitude_s, r=(N/1000)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  scale_fill_manual(values=colz)


# 
# # copie of latitude (x.) and longitude (y.)
pies$x2 <- pies$longitude_e
pies$y2 <- pies$latitude_s
# 
# #manually move the pie chart in case there is an overlap (change x and y)


pies$x2[pies$bat_species== "Eidolon dupreanum"] <- pies$longitude_e[pies$bat_species== "Eidolon dupreanum"] + 4
pies$y2[pies$bat_species== "Eidolon dupreanum"] <- pies$latitude_s[pies$bat_species== "Eidolon dupreanum"] - 1

head(pies)



#plot pie chart 
#loko<-c("Rousettus madagascariensis"="#B200ED","Eidolon dupreanum"="#7FFF00","Pteropus rufus"="#0000FF")

#this is Fig1A
Fig1a <- p2b+
  annotate("segment", x=pies$longitude_e, xend=pies$x2,y=pies$latitude_s,yend=pies$y2,size=.7)+ # put the lines
  geom_scatterpie(aes(x=x2, y=y2, r=(log10(N)/1.2)), 
                  data = pies, cols="plot_class", long_format=TRUE) +
  theme_bw() +theme(panel.grid = element_blank(),
                    plot.title = element_text(color="black", size=12, face="bold"),
                    plot.margin = unit(c(0.7,.5,.3,.1),"cm"),
                    axis.title.x = element_text(color="black", size=12),
                    axis.title.y = element_text(color="black", size=12),
                    legend.position=c(.85,.85),
                    legend.margin = margin(),
                    legend.title=element_blank(),
                    legend.text = element_text(size = 7.5)) +
  scale_fill_manual(values=colz) +
  coord_sf(xlim=c(41,54))#+
  #geom_scatterpie_legend(log10(c(10,100)/1.2),
   #                      x=54.5, y=-23.5, 
    #                     n=2,
     #                    labeller = function(x) paste(10^(x)*1.2,"indiv"))

#print(p4)


ggsave(file = paste0(homewd, "final-figures/tmp.pdf"),
       plot=Fig1a,
       units="mm",  
       width=60, 
       height=65, 
       scale=3, 
       dpi=300)







#and this is Fig 1B

#get into date form
dat$collection_date <- as.Date(dat$collection_date,format = "%m/%d/%y")

#check sites
unique(dat$roost_site) #all sites are in the moramanga area and can be treated as one

#now there are some urine and some fecal samples
#one of the urine samples is the same individual as one of the fecal samples
#so go ahead and remove the fecal sample (negative) from consideration



#get the date of the first day of every week
dat$epiwk <- cut(dat$collection_date, "week")
dat$epiwk <- as.Date(as.character(dat$epiwk))


names(dat)[names(dat)=="bat_species"] <- "species"

#and make sure it is only 1 of the same sample type per each date

dat.list <- dlply(dat, .(sampleid))

#for those with multiple samples, slim to one
out = c(unlist(lapply(dat.list, nrow)))
out[out>1] #none


#should be two urine samples
length(dat$sample_type[dat$sample_type=="urine"])#106

#okay to go.

#summarize into prevalence by species and epiwk
dat.sum <- ddply(dat, .(species, epiwk), summarise, N=length(henipavirus), pos=sum(henipavirus))

#get negatives and prevalence
dat.sum$neg= dat.sum$N-dat.sum$pos
dat.sum$prevalence <- dat.sum$pos/dat.sum$N

#and confidence intervals on the prevalence
CIs <- mapply(FUN=prop.test, x=as.list(dat.sum$pos), n=as.list(dat.sum$N), MoreArgs = list(alternative = "two.sided", conf.level = .95, correct=F), SIMPLIFY = F)

#and extract the upper and lower CIs
get.CIs <- function(df){
  lci = df$conf.int[1]
  uci = df$conf.int[2]
  out.df<- cbind.data.frame(lci=lci, uci=uci)
  return(out.df)
}

CIs <- lapply(CIs, get.CIs)

dat.sum$lci <- c(unlist(sapply(CIs, '[',1)))
dat.sum$uci <- c(unlist(sapply(CIs, '[',2)))

#simplify= the name of "host_genus_species" 
names(dat.sum)[names(dat.sum)=="host_genus_species"] <- "species"

#and plot
#here's a vector assigning colors to each species
colz = c("Eidolon dupreanum"="steelblue1")#, "Pteropus rufus" = "violetred", "Rousettus madagascariensis" = "seagreen" )
names(dat.sum)[names(dat.sum)=="age_class"] <- "age"
# dat.sum$age[dat.sum$age=="A"] <- "adult"
# dat.sum$age[dat.sum$age=="J"] <- "juvenile"

#shapez = c("juvenile" = 17, "adult" = 16)

p1 <- ggplot(data=dat.sum) + #here is the dataset
  geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci, color=species), size=.1) + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=epiwk, y= prevalence, color=species, size=N)) + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("henipavirus prevalence") + #change the name of the y-axis
  scale_color_manual(values=colz) + #assign the colors manually using the vector above
  #scale_shape_manual(values=shapez) +
  #facet_grid(age~.) +
  theme_bw() + #some style features
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic"),
        strip.background = element_rect(fill="white"),
        axis.title.x = element_blank()) 


#p1

#most information seems to be captured in plot 1, so lets keep that

#or try another plot that includes the winter season for Moramanga instead
seas.dat1 = cbind.data.frame(x=as.Date(c("2018-06-01","2018-09-01")), ymin=c(-1,-1), ymax=c(2,2), label="dry season")
seas.dat2 = cbind.data.frame(x=as.Date(c("2018-02-01","2018-06-01")), ymin=c(-1,-1), ymax=c(2,2), label="late-stage juveniles")

#jitter dates manually
# dat.sum$epiwk_jitter <- dat.sum$epiwk
# #dat.sum$epiwk_jitter[dat.sum$species=="Rousettus madagascariensis"] <- dat.sum$epiwk_jitter[dat.sum$species=="Rousettus madagascariensis"] + 1
# dat.sum$epiwk_jitter[dat.sum$species=="Eidolon dupreanum"] <- dat.sum$epiwk_jitter[dat.sum$species=="Eidolon dupreanum"] + 2

text.dat = cbind.data.frame(date=c(as.Date("2018-07-15")),y=c(.97), label=c( "dry season"))

Fig1b <-  ggplot(data=dat.sum) +
  geom_ribbon(data=seas.dat1, aes(x=x, ymin=ymin, ymax=ymax), fill="cornflowerblue", alpha=.3) +
  #geom_ribbon(data=seas.dat2, aes(x=x, ymin=ymin, ymax=ymax), fill="lightgoldenrod1", alpha=.3) +
  geom_errorbar(aes(x=epiwk, ymin=lci, ymax=uci), size=.2, color="mediumseagreen") + #here we plot the uci and lci for prevalence, with lines colored by species
  geom_point(aes(x=epiwk, y= prevalence, color=species, size=N), color="mediumseagreen") + #here we plot the mean prevalence, with dot colored by species and sized based on sample size per date
  ylab("henipavirus prevalence") + #change the name of the y-axis
  geom_text(data=text.dat, aes(x=date,y=y, label=label), size=3) +
  #geom_text(data=text.dat, aes(x=as.Date("2018-04-24"),y=.95, label="late-stage juveniles"), size=3) +
  #scale_color_manual(values=colz) + #assign the colors manually using the vector above
  #scale_fill_manual(values=colz) + #assign the colors manually using the vector above
  #scale_shape_manual(values=shapez) +#
  theme_bw() + #some style features
  theme(panel.grid = element_blank(), legend.text = element_text(face="italic", size = 8),
        legend.position = c(.1,.8),
        strip.background = element_rect(fill="white"),
        legend.box = "horizontal",
        legend.background = element_rect(fill="white"),
        #legend.title = element_blank(),
        legend.spacing.y = unit(.05, "cm"),
        plot.margin = unit(c(.8,.2,.7,.3),"cm"),
        axis.title.x = element_blank()) +  #more style features
        coord_cartesian(xlim=as.Date(c("2018-01-01", "2019-05-01")), ylim=c(0,1))
#Fig1b



Fig1all <- cowplot::plot_grid(Fig1a, Fig1b, nrow=1, ncol=2, labels = c("(A)", "(B)"), label_x = c(.1,-.01), label_y = .99)


ggsave(file = paste0(homewd, "final-figures/Fig1.pdf"),
       plot=Fig1all,
       units="mm",  
       width=150, 
       height=65, 
       scale=3, 
       dpi=300)




