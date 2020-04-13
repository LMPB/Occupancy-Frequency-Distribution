#Script for the map of world distribution of core bacteria
  #Library upload
  library(ggplot2)
  library(vegan)
  library(tidyr)
  library(dplyr)
  #Map
  library(rgdal)
  library(raster)
  library(broom)
  library(RColorBrewer)
  library(dplyr)
  library(maps) #mapas simples, eixos, escala, cidades 
  library(mapdata) #base de dados WorldHires e rios
  library(rworldmap) #outra base de dados de mapas do mundo
  library(maptools) #Ler ESRI shapefiles 
  library(mapproj) #Projeções e grids
  library(ggmap) #Gmaps, OSM + mapas baseados em ggplot2
  library(rgdal)
  library(remotes)
  library(grid)
  #grid all the maps
  library (cowplot)
  library(gridExtra)
  library(ggpubr)

#Data organizing
#In this step we had two datasets: asv_global and asv_global_count
#asv_global have 7 collumns:
  #asv: sequence name
  #location: lake or country name
  #region: continents
  #lat: latitude degree
  #long: longitude degree
  #isolation_source: type of environment (e.g. lake, reservoir,river)
  #ref_group: taxonomic classification (e.g. hgcI clade, Polynucleobacter)
  
#asv_global_count have 4 columns:
  #location: lake or country name	
  #count_location: times a sequence was found in each location	
  #lat: latitude degree
  #long: longitude degree
  
setwd("")

core2 <- read.table("asv_global.txt", header=T,sep="\t", dec=".") 
count <- read.table("asv_global_count.txt", header=T, sep="\t", dec=".")
latlong2= core2
#core2 <- core2[2:5]
#core2 <- core2[-3]
core.count = count (core2, vars = isolation_source)
core.count.region = count (core2, vars = region)

core.count$vars <- factor(core.count$vars, levels = c("lake","reservoir","river","drinking water","floodplain","groundwater","estuary","other"))
core.count.region$vars <- factor(core.count.region$vars, levels = c("Africa", "Asia","Oceania","Europe","South America","Central America","North America"))
latlong2$ref_group <- factor(latlong2$ref_group, levels = c("hgcI_clade", "Polynucleobacter" ))
#count$ref_group <- factor(count$ref_group, levels = c("hgcI_clade", "Polynucleobacter" ))
#count$ref_group <- factor(count$ref_group)


#-------------------
#General global distribution
#Creating the pie graph -----------------------------
system =
  ggplot(core.count, aes(x=" ", y=n, fill=vars))+ 
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+ 
  geom_text(aes(label = paste0(round(n*1))), size=3, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values=c("seagreen3", "steelblue2","olivedrab2", "#33658A", "#2F4858", "#F6AE2D", "#F26419", "#999999")) + 
  theme(axis.line = element_blank(), axis.ticks = element_blank())+
  labs(fill="System Type", x=NULL, y=NULL, title=" ")+
  theme_bw()+
  theme(legend.position = "")
system

#creating the base map
world_map <- map_data("world")

p <- ggplot() + 
  coord_fixed() +
  xlab("") + 
  ylab("")

base_world_messy <- p + 
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
               colour="grey50", fill="grey50")

#Ploting the data into the base map
map =
  base_world_messy + 
  geom_point(data=count,aes(x=long, y=lat), size = 3, alpha=I(0.7)) +
  scale_color_manual (values=c("red","grey20","purple"))+ 
  labs( x="Longitude (degrees)", y="Latitude (degrees)")+
  theme_bw()+
  theme(legend.position = "")
map

#set the pie graph inside the world map
#manipulating the xmin, xmax,ymin and ymax values we can determine the location were the pie graph will be set
map.compl = map + 
  annotation_custom(ggplotGrob(system), xmin = 50, xmax = 110, ymin = -85, ymax = 0)

map.compl

#-------------------
#Global distribution of each ASV

#ASV_00001
core.ASV_00001 = filter(core2,asv=="ASV_00001") 
core.count.ASV_00001 = count (core.ASV_00001, vars = isolation_source)
core.count.ASV_00001 = transform(core.count.ASV_00001, vars=factor(vars,levels = c("lake", "reservoir", "river","floodplain","estuary","other")))

latlong.ASV_00001 = filter(latlong2,asv=="ASV_00001")

#asv 1 pie graph
system.ASV_00001 =
  ggplot(core.count.ASV_00001, aes(x=" ", y=n, fill=vars))+ 
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+ 
  geom_text(aes(label = paste0(round(n*1))), size=3, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values=c("seagreen3","steelblue2","olivedrab2","#2F4858","#F26419",  "#999999")) + 
  theme(axis.line = element_blank(), axis.ticks = element_blank())+
  labs(fill="System Type", x=NULL, y=NULL, title=" ")+
  theme_bw()+
  theme(legend.position = "")


#creating the base map
world_map <- map_data("world")

p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                     colour="grey50", fill="grey50")

#Ploting the data into the base map
map.ASV_00001 =
  base_world_messy + 
  geom_point(data=latlong.ASV_00001,aes(x=long, y=lat, colour=ref_group), alpha=I(0.7), size = 2) +
  scale_color_manual (values=c("red"))+ 
  labs( x="Longitude (degrees)", y="Latitude (degrees)")+
  theme_bw()+
  ggtitle("ASV_00001")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "")

#set the pie graph inside the world map
map.compl.ASV_00001 = map.ASV_00001 + 
  annotation_custom(ggplotGrob(system.ASV_00001), xmin = -210, xmax = -90, ymin = -70, ymax = 15)

map.compl.ASV_00001

#ASV_00002
core.ASV_00002 = filter(core2,asv=="ASV_00002") 
core.count.ASV_00002 = count (core.ASV_00002, vars = isolation_source)
core.count.ASV_00002 = transform(core.count.ASV_00002, vars=factor(vars,levels = c("lake","reservoir","river","drinking water","floodplain","estuary","other")))

latlong.ASV_00002 = filter(latlong2,asv=="ASV_00002")

#asv 2 pie graph
system.otu_2 =
  ggplot(core.count.ASV_00002, aes(x=" ", y=n, fill=vars))+ 
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+ 
  geom_text(aes(label = paste0(round(n*1))), size=3, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values=c("seagreen3", "steelblue2","olivedrab2", "#33658A", "#2F4858", "#F26419", "#999999")) + 
  theme(axis.line = element_blank(), axis.ticks = element_blank())+
  labs(fill="System Type", x=NULL, y=NULL, title=" ")+
  theme_bw()+
  theme(legend.position = "")

#creating the base map
world_map <- map_data("world")

p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                     colour="grey50", fill="grey50")


#Ploting the data into the base map
map.ASV_00002 =
  base_world_messy + 
  geom_point(data=latlong.ASV_00002,aes(x=long, y=lat, colour=ref_group), alpha=I(0.7), size = 2) +
  scale_color_manual (values=c("red"))+ 
  labs( x="Longitude (degrees)", y="Latitude (degrees)")+
  theme_bw()+
  ggtitle("ASV_00002")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "")

#set the pie graph inside the world map
map.compl.ASV_00002 = map.ASV_00002 + 
  annotation_custom(ggplotGrob(system.ASV_00002), xmin = -210, xmax = -90, ymin = -70, ymax = 15)

map.compl.ASV_00002

#ASV_00003
core.ASV_00003 = filter(core2,asv=="ASV_00003") 
core.count.ASV_00003 = count (core.ASV_00003, vars = isolation_source)
core.count.ASV_00003 = transform(core.count.ASV_00003, vars=factor(vars,levels = c("lake","reservoir","river","drinking water","floodplain","groundwater","other")))

latlong.ASV_00003 = filter(latlong2,asv=="ASV_00003")

#asv 3 pie graph
system.ASV_00003 =
  ggplot(core.count.ASV_00003, aes(x=" ", y=n, fill=vars))+ 
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+ 
  geom_text(aes(label = paste0(round(n*1))), size=3, position = position_stack(vjust = 0.5))+
  scale_fill_manual(values=c("seagreen3", "steelblue2","olivedrab2", "#33658A", "#2F4858", "#F6AE2D", "#999999")) + 
  theme(axis.line = element_blank(), axis.ticks = element_blank())+
  labs(fill="System Type", x=NULL, y=NULL, title=" ")+
  theme_bw()+
  theme(legend.position = "")

#creating the base map
world_map <- map_data("world")

p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                     colour="grey50", fill="grey50")

#Ploting the data into the base map
map.ASV_00003 =
  base_world_messy + 
  geom_point(data=latlong.ASV_00003,aes(x=long, y=lat, colour=ref_group), alpha=I(0.7), size = 2) +
  scale_color_manual (values=c("darkgreen"))+ 
  labs( x="Longitude (degrees)", y="Latitude (degrees)")+
  theme_bw()+
  ggtitle("ASV_00003")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "")

#set the pie graph inside the world map
map.compl.ASV_00003 = map.ASV_00003 + 
  annotation_custom(ggplotGrob(system.ASV_00003), xmin = -210, xmax = -90, ymin = -70, ymax = 15)

map.compl.ASV_00003

#Concatenate all the graphs provided above in one figure
library(gridExtra)
library(ggpubr)
library(cowplot)

#Arrange the figures
x = arrangeGrob(map.compl,
                map.ASV_00001,map.ASV_00002,
                map.ASV_00003,
                ncol = 2, nrow = 3,
                layout_matrix = rbind(c(1,1), c(2,3), c(4,5)))

#set 
as_ggplot(x) +
  draw_plot_label(label = c ("A","B","C","D","E"), size = 15,
                  x = c(0, 0, 0.5,0,0.5), y = c(1, 0.60, 0.60,0.25,0.25))

#######################################
#Sampling location map

#Packages
library(maps)  
library(mapdata) 
library(rworldmap) 
library(maptools)  
library(mapproj) 
library(ggmap) 
library(rgdal)
library(remotes)
library(grid)

#Dowload shapefiles
maptools::readShapePoly
rgdal::readOGR

#to create a sampling map, download *.shp file from institutional websites
#Getting the Latin American map from mapdata package
la.map = borders("worldHires", 
                           regions = c("Brazil", "Uruguay", "Argentina", "French Guiana", "Suriname", "Colombia", 
                                                     "Venezuela", "Bolivia", "Ecuador", "Chile", "Paraguay", "Peru", "Guyana", 
                                                     "Panama", "Costa Rica", "Nicaragua", "Honduras", "El Salvador", "Belize", 
                                                     "Guatemala", "Mexico", "Trinidad and Tobago", "Caribe", "Puerto Rico", 
                                                     "Dominican Republic", "Haiti", "Jamaica", "Cuba", "Bahamas", "Antiles", 
                                                     "Dominica", "Saba"), fill = "grey50", colour = "grey50")

#getting a region from the country shp archive
br = readOGR("my/documents/BRASIL.shp")
proj4string(br)
names(br)
head(br$Estado)
sp = br[br@data$Estado =="São Paulo",] 
sp.f = fortify(sp)

#setting the geographical location of each sampling site
#geo1 have 3 collumns:
  #otu: sequence names
  #lat: latitude degree
  #long: longitude degree
geo <- read.table("my/documents/geo1.txt", header=TRUE, sep="\t", 
                  na.strings="NA", dec=".", strip.white=TRUE)


# Scale Bar Coordinates Function ------------------------------------------
#downloaded from https://stackoverflow.com/questions/17151024/adding-scale-bar-to-ggplot-map

# Return a list whose elements are :
# 	- rectangle : a data.frame containing the coordinates to draw the first rectangle ;
# 	- rectangle2 : a data.frame containing the coordinates to draw the second rectangle ;
# 	- legend : a data.frame containing the coordinates of the legend texts, and the texts as well.
#
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distance_lon : length of each rectangle ;
# distance_lat : width of each rectangle ;
# distance_legend : distance between rectangles and legend texts ;
# dist_units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
create_scale_bar <- function(lon,lat,distance_lon,distance_lat,distance_legend, dist_units = "km"){
  # First rectangle
  bottom_right <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon, dist.units = dist_units, model = "WGS84")
  
  topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_lat, dist.units = dist_units, model = "WGS84")
  rectangle <- cbind(lon=c(lon, lon, bottom_right[1,"long"], bottom_right[1,"long"], lon),
                     lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
  rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
  
  # Second rectangle t right of the first rectangle
  bottom_right2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon*2, dist.units = dist_units, model = "WGS84")
  rectangle2 <- cbind(lon = c(bottom_right[1,"long"], bottom_right[1,"long"], bottom_right2[1,"long"], bottom_right2[1,"long"], bottom_right[1,"long"]),
                      lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
  rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
  
  # Now let's deal with the text
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_legend, dist.units = dist_units, model = "WGS84")
  on_top2 <- on_top3 <- on_top
  on_top2[1,"long"] <- bottom_right[1,"long"]
  on_top3[1,"long"] <- bottom_right2[1,"long"]
  
  legend <- rbind(on_top, on_top2, on_top3)
  legend <- data.frame(cbind(legend, text = c(0, distance_lon, distance_lon*2)), stringsAsFactors = FALSE, row.names = NULL)
  return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}


# North Arrow Function ----------------------------------------------------
# Returns a list containing :
#	- res : coordinates to draw an arrow ;
#	- coordinates of the middle of the arrow (where the "N" will be plotted).
#
# Arguments : #
#-------------#
# scale_bar : result of create_scale_bar() ;
# length : desired length of the arrow ;
# distance : distance between legend rectangles and the bottom of the arrow ;
# dist_units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
create_orientation_arrow <- function(scale_bar, length, distance = 1, dist_units = "km"){
  lon <- scale_bar$rectangle2[1,1]
  lat <- scale_bar$rectangle2[1,2]
  
  # Bottom point of the arrow
  beg_point <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance, dist.units = dist_units, model = "WGS84")
  lon <- beg_point[1,"long"]
  lat <- beg_point[1,"lat"]
  
  # Let us create the endpoint
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = length, dist.units = dist_units, model = "WGS84")
  
  left_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 225, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  right_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 135, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  res <- rbind(
    cbind(x = lon, y = lat, xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = left_arrow[1,"long"], y = left_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = right_arrow[1,"long"], y = right_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]))
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  # Coordinates from which "N" will be plotted
  coords_n <- cbind(x = lon, y = (lat + on_top[1,"lat"])/2)
  
  return(list(res = res, coords_n = coords_n))
}


# Draw Scale Bar and North Arrow  -----------------------------------------

#
# Result #
#--------#
# This function enables to draw a scale bar on a ggplot object, and optionally an orientation arrow #
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distance_lon : length of each rectangle ;
# distance_lat : width of each rectangle ;
# distance_legend : distance between rectangles and legend texts ;
# dist_units : units of distance "km" (kilometers) (by default), "nm" (nautical miles), "mi" (statute miles) ;
# rec_fill, rec2_fill : filling colour of the rectangles (default to white, and black, resp.);
# rec_colour, rec2_colour : colour of the rectangles (default to black for both);
# legend_colour : legend colour (default to black);
# legend_size : legend size (default to 3);
# orientation : (boolean) if TRUE (default), adds an orientation arrow to the plot ;
# arrow_length : length of the arrow (default to 500 km) ;
# arrow_distance : distance between the scale bar and the bottom of the arrow (default to 300 km) ;
# arrow_north_size : size of the "N" letter (default to 6).
scale_bar <- function(lon, lat, distance_lon, distance_lat, distance_legend, dist_unit = "km", rec_fill = "white", rec_colour = "black", rec2_fill = "black", rec2_colour = "black", legend_colour = "black", legend_size = 3, orientation = TRUE, arrow_length = 500, arrow_distance = 300, arrow_north_size = 6){
  the_scale_bar <- create_scale_bar(lon = lon, lat = lat, distance_lon = distance_lon, distance_lat = distance_lat, distance_legend = distance_legend, dist_unit = dist_unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = the_scale_bar$rectangle, aes(x = lon, y = lat), fill = rec_fill, colour = rec_colour)
  # Second rectangle
  rectangle2 <- geom_polygon(data = the_scale_bar$rectangle2, aes(x = lon, y = lat), fill = rec2_fill, colour = rec2_colour)
  # Legend
  scale_bar_legend <- annotate("text", label = paste(the_scale_bar$legend[,"text"], dist_unit, sep=""), x = the_scale_bar$legend[,"long"], y = the_scale_bar$legend[,"lat"], size = legend_size, colour = legend_colour)
  res <- list(rectangle1, rectangle2, scale_bar_legend)
  
  if(orientation){# Add an arrow pointing North
    coords_arrow <- create_orientation_arrow(scale_bar = the_scale_bar, length = arrow_length, distance = arrow_distance, dist_unit = dist_unit)
    arrow <- list(geom_segment(data = coords_arrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coords_arrow$coords_n[1,"x"], y = coords_arrow$coords_n[1,"y"], size = arrow_north_size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res)
}


# Draw the maps ------------------------------------------------------------
map.la = 
  
  ggplot() +
  #geom_point(data = geo, aes(x = long, y = lat, group = region, color = region), size = 1)+
  la.map+
  #geom_polygon(data = bra.map, aes(x = long, y = lat, group = group)) +
  coord_equal()+ # map_data cria um dataframe que deve ser acrescentado como geom_polygon()
  scale_bar(lon = -74, lat = 30, distance_lon = 2000, distance_lat = 200, distance_legend = 550, dist_unit = "km", orientation = FALSE)+
  geom_polygon(data = sp1, aes(x = long, y = lat, group = group),fill = "gray40", color = "black")+
  #geom_polygon(data = par.f, aes(x = long, y = lat, group = group, fill = "lightblue"), alpha = 0.3)+
  #geom_point(data = geo, aes(x = long, y = lat, group = region, color = region))+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = -23.5, linetype = "dashed")+
  #scale_color_manual(name = "Sub Basin", labels=c("Grande","Paranapanema","Serra do mar","Tietê"), values = c("orange","red","blue","green"))+
  scale_fill_manual (name = "", labels = ("Paraná Basin"), values = ("lightblue"))+
  scale_color_manual(name = "Sub Basin", 
                     labels=c("Grande","Paranapanema","Serra do mar","Tietê"), 
                     values = c("orange","red","blue","green"))+
  theme_bw()+
  ylab("Latitude")+
  xlab("Longitude")+
  geom_text(show.legend = FALSE)+
  theme(legend.position = "none")

map.la


map.sp =
  
  ggplot()+
  geom_polygon(data = sp.f, aes(x = long, y = lat, group = group), fill = "gray40", color = "black") +
  #geom_polygon(data = par.f, aes(x = long, y = lat, group = group, fill = "lightblue"), alpha = 0.3)+
  coord_equal()+ # map_data cria um dataframe que deve ser acrescentado como geom_polygon()
  scale_bar(lon = -53, lat = -25, distance_lon = 100, distance_lat = 30, distance_legend = 50, dist_unit = "km", orientation = FALSE)+
  geom_point(data = geo, aes(x = long, y = lat, group = region), color = "black", size = 3)+
  geom_hline(yintercept = -23.5, linetype = "dashed")+
  scale_color_manual(name = "region", labels=c("Grande","Paranapanema","Serra do mar","Tietê"), values = c("green", "orange","red","blue","green"))+
  geom_text(show.legend = FALSE)+
  theme_bw()+
  ylab("Latitude")+
  xlab("Longitude")
  #theme(legend.position = "none")

map.sp

#Grob the Latin American graph inside the regional map
map.compl = map.la + 
  annotation_custom(ggplotGrob(map.sp), xmin = -40, xmax = -20, ymin = -38, ymax = -18)

map.compl


