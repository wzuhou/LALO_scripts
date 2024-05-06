{library(tidyverse)
library(ggplot2)
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
setwd('M:/ROSLIN/LALO_genome/Map/')
}
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

##############
#add geo distribution
#Birdlife 2017 data
library(rgdal)
# The input file geodatabase-ALL birds
gdb_file = paste0( "./Bird_Map/BOTW_2022_1/BOTW.gdb")
# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list = ogrListLayers(gdb_file)
print(fc_list)
botw <-st_read(dsn=gdb_file, layer="All_Species")
names(botw)

#LALO
LALO <- subset(botw,sisid=='22721033') #22721088 -WCS
st_write(LALO, "./Bird_Map/BOTW_2022_1/shapefile/Lalo.shp", driver="ESRI Shapefile")
class(LALO)
df <- LALO[,'seasonal']

df_union_cast <- st_cast(df, "POLYGON")
#this is the right one!
df_union_cast$seasonal<-as.factor(df_union_cast$seasonal)


#################Add points-station
station <- tibble::tribble( 
  ~city,           ~lat,     ~lon, ~colour, ~ Shape
)

p1 <- ggplot(data = world) +
 #geom_map(data = world, map = world,aes(long, lat, map_id = region),color = "gray40", fill = "gray80", size = 0.1) +
  geom_sf(fill= 'ghostwhite',size=0.25,colour='gray40' ) +
  geom_sf(aes(fill = seasonal),data=df_union_cast, alpha=0.7,colour= 'gray70',size=0.1)+
  geom_point(data = station[1,], mapping = aes(x = lon, y = lat),size=1.5, shape=21,fill ='#972D15',alpha=0.9,color='#972D15') +
  #scale_x_continuous(breaks = seq(-60,-160,by=-20))+
  scale_y_continuous(breaks = seq(20,100,20))+
  xlab("Longitude") + ylab("Latitude") + 
  scale_fill_manual(values = c('#02401B','#D8B70A','#A2A475'),labels = c("Breeding", "Non-breeding",'Passage'))+
  #scale_colour_manual(values = c('#DD8D29','#46ACC8','#B40F20'))+
  #annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-170, -50.12), ylim = c(10, 85))+
  guides(colour = guide_legend(title ="Seasonal",keyheight = 0.8,label.vjust=1,ncol = 1))+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = linewidth(fill = NULL),
        legend.key.size = unit(0.3, 'cm'),
        legend.background = element_rect(size=0.5)
        )

p1

ggsave('Map_North_America_distribution.tiff',p1,width = 5,height =7 , dpi = 350)#width = 5,height =7 


###Zoomed figure
p2 <- ggplot(data = world) +
    geom_sf(fill= 'ghostwhite',size=0.3 ) +
  xlab(NULL) + ylab(NULL) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering) +

  geom_point(data = station, mapping = aes(x = lon, y = lat),size=2.5, shape=station$Shape,fill =station$colour,alpha=0.8,color=station$colour) +
  coord_sf(xlim = c(-160, -140), ylim = c(65, 73))+
  #coord_sf(xlim = c(-170, -50.12), ylim = c(10, 85))+
  scale_y_continuous(breaks = seq(65,80,5))+
  scale_x_continuous(breaks = seq(-140,-160,by=-20))+
  annotation_scale(location = "tl", width_hint = 0.2,pad_x = unit(0.05, "in"), pad_y = unit(0.6, "in")) +
  annotate(geom = "text", x = -144, y = 71, label = "Arctic Ocean", fontface = "bold", color = "grey22", size = 6) +
  annotate(geom = "text", x = -144, y = 70.5, label = "Beaufort Sea", fontface = "italic", color = "grey22", size = 4) +
  annotate(geom = "text", x = -150, y = 65, label = "Alaska", fontface = "italic", color = "darkblue", size = 4) +
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = NULL)
  )

p2
ggsave('Map_Sampling.tiff',p2, width = 5,height =5.2 ,dpi = 300)

###END###
