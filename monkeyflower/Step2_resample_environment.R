####Script to resample the environmental data to 9m,27m,81m
####Author: Ian Breckheimer
####2 November 2013

library(raster)
library(rgdal)

####Brings in the environmental variables####
setwd("/Users/ian/GIS/")
env3 <- brick("env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)")
names(env3) <- env_names

##Resamples environmental data.
env9 <- aggregate(env3,fun=mean,na.rm=T,fact=3,
                  filename="env_pred_9m.img",overwrite=T)
names(env9) <- env_names
env27 <- aggregate(env3,fun=mean,na.rm=T,fact=9,
                   filename="env_pred_27m.img",overwrite=T)
names(env27) <- env_names
env81 <- aggregate(env3,fun=mean,na.rm=T,fact=27,
                   filename="env_pred_81m.img",overwrite=T)
names(env81) <- env_names

# ##Examines the resampled data.
# disp <- stack(env24[[c(4,7,3,10)]])
# disp_titles <- env_titles[c(4,7,3,10)]
# par(mar = c(0.1, 0.1, 1, 0.1), mfrow = c(2,2))
# pals <- list(terrain.colors(255),
#              rainbow(255,start=0.3,end=0.9),
#              rev(heat.colors(255)),
#              rainbow(255,start=0.05,end=0.7))
# for (i in 1:(nlayers(disp)-1)) {
#   plot(disp[[i]],axes=FALSE,legend=T,asp=1,box=F,col=pals[[i]])
#   title(main=disp_titles[i],line=-0.5)
# }
# plot(disp[[4]],axes=FALSE,legend=T,asp=1,box=F,col=pals[[4]])
# title(main=disp_titles[[4]],line=-0.5)
# axis(1,line=-2,at=c(592000,602000),labels=c("0","10km"))