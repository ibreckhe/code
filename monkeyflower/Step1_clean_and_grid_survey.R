##Imports mimulus survey data, merges high and low-density surveys, and grids the data.
##Ian Breckheimer
##25 October 2013

##Sets up the workspace
library(data.table)
library(plyr)
library(maptools)
library(raster)
library(rgdal)

##Brings in the 3m environmental data to serve as a raster template.
setwd("/Users/ian/GIS/")
env3 <- brick("env_pred_3m.img",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")
env_titles <- c("Topographic Roughness","Canopy Ht. (m)","Canopy Cover (prop.)","Elevation (m)","Insolation Time (hrs)",
                "Slope (deg)","Solar Rad (Wh/sq.m)","Topo Wetness","Max Topo Wetness","Snow Melt","NDVI",
                "Lidar Intensity","Pond-Stream Area (sq m)","Stream Distance (m)")
names(env3) <- env_names

####Merges high and low density survey data####
setwd("/Users/ian/Dropbox/Research/monkeyflower/data/survey_data/")

##Reads the data into R as dataframes.
sites <- read.table("Survey_data_highdensity_2013.csv",sep=",",header=T)
patches <- read.table("Survey_data_high_density_2013_patches.csv",sep=",",header=T)
sites$UID <- sites$meta.instanceID
patches$UID <- patches$parent_uid

##Converts them to data.tables.
sites_dt <- data.table(sites,key="UID")
patches_dt <- data.table(patches,key="UID")

##Merges the two data frames based on the uuid field (too slow)
all_dt <- merge(sites_dt,patches_dt,all.y=TRUE)

##Sums all of the patch attributes in every site using plyr.
all_sum <- ddply(all_dt, c("UID"), function(df) return(c(Stem_sum=sum(df$Stem_num,na.rm=T),
                                                         Area_sum_cm=sum(df$Patch_extent_cm,na.rm=T),
                                                         Flower_sum=sum(df$flower_num,na.rm=T),
                                                         Num_patches=length(df$Patch_extent_cm))))
all <- cbind(sites_dt,all_sum)
rownames(all) <- all$UID

##Converts the data to a SpatialLinesDataFrame.
line_fun <- function(df){Lines(Line(cbind(c(df$location_start.Longitude,df$location_end.Longitude),
                                          c(df$location_start.Latitude,df$location_end.Latitude))),ID=df$UID)}
all_line <- dlply(all,"UID",line_fun)
all_spatial <- SpatialLines(all_line, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
all_spdf <- SpatialLinesDataFrame(all_spatial,all)

##Writes the data frame to a shapefile.
writeSpatialShape(all_spdf, "survey_highdensity_lines", factor2char = TRUE, max_nchar=254)

####Grids the cleaned survey data####
setwd("/Users/ian/Dropbox/Research/monkeyflower/data/survey_data/")

##Imports the cleaned data from shapefile format.
g_lowd <- readShapePoints("cleaned/guttatus_lowdensity_points",proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
g_highd <- readShapePoints("cleaned/guttatus_highdensity_points",proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
t_lowd <- readShapePoints("cleaned/tilingii_lowdensity_points",proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))
t_highd <- readShapePoints("cleaned/tilingii_highdensity_points",proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))

track <- readShapeLines("cleaned/survey_tracks_2013_singlepart",proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"))

##Transforms to UTM coordinates.
g_lowd_utm10 <- spTransform(g_lowd,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
g_highd_utm10 <- spTransform(g_highd,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
t_lowd_utm10 <- spTransform(t_lowd,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
t_highd_utm10 <- spTransform(t_highd,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
track_utm10 <- spTransform(track,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

##Rasterizes the point count field.
setwd("/Users/ian/GIS/tmp/")
g_lowd_ras <- rasterize(g_lowd_utm10,env3,field="Stemnum",fun=sum,filename="g_lowd_3m.tif",overwrite=T)
g_highd_ras <- rasterize(g_highd_utm10,env3,field="Stem_sum",fun=sum,filename="g_highd_3m.tif",overwrite=T)
t_lowd_ras <- rasterize(t_lowd_utm10,env3,field="Stemnum",fun=sum,filename="t_lowd_3m.tif",overwrite=T)
t_highd_ras <- rasterize(t_highd_utm10,env3,field="Stem_sum",fun=sum,filename="t_highd_3m.tif",overwrite=T)

##Creates a layer representing all cells sampled.
allfun <- function(x){if (sum(x,na.rm=T) > 0) return(1) else return(NA)}
all_samp <- calc(stack(g_lowd_ras,g_highd_ras,t_lowd_ras,t_highd_ras,nod_ras),
                 fun=allfun,filename="all_sampled_3m.tif",overwrite=T)

##Adds low and high-density surveys together.
gut_data <- stack(g_lowd_ras,g_highd_ras,all_samp)
til_data <- stack(t_lowd_ras,t_highd_ras,all_samp)
gut_tot <- calc(gut_data,fun=function(x){sum(x,na.rm=T)-1},filename="gut_tot_3m.tif",overwrite=T)
til_tot <- calc(til_data,fun=function(x){sum(x,na.rm=T)-1},filename="til_tot_3m.tif",overwrite=T)

##Rasterizes the patch extent field.
g_lowd_ext_ras <- rasterize(g_lowd_utm10,env3,field="Patchexten",fun=sum,filename="g_lowd_ext_3m.tif",overwrite=T)
g_highd_ext_ras <- rasterize(g_highd_utm10,env3,field="Area_sum_c",fun=sum,filename="g_highd_ext_3m.tif",overwrite=T)
t_lowd_ext_ras <- rasterize(t_lowd_utm10,env3,field="Patchexten",fun=sum,filename="t_lowd_ext_3m.tif",overwrite=T)
t_highd_ext_ras <- rasterize(t_highd_utm10,env3,field="Area_sum_c",fun=sum,filename="t_highd_ext_3m.tif",overwrite=T)

##Adds low and high-density surveys together.
gut_ext_data <- stack(g_lowd_ext_ras,g_highd_ext_ras,all_samp)
til_ext_data <- stack(t_lowd_ext_ras,t_highd_ext_ras,all_samp)
gut_ext_tot <- calc(gut_ext_data,fun=function(x){sum(x,na.rm=T)-1},filename="gut_ext_3m.tif",overwrite=T)
til_ext_tot <- calc(til_ext_data,fun=function(x){sum(x,na.rm=T)-1},filename="til_ext_3m.tif",overwrite=T)

##Puts them all in a raster brick.
mim <-  stack(gut_tot,til_tot,gut_ext_tot,til_ext_tot)
mim_brick <- brick(mim,filename="mim_survey_2013_multiband.tif",overwrite=T)
names(mim_brick) <- c("gut_stems","til_stems","gut_ext_cm","til_ext_cm")
mim_brick_zero <- reclassify(mim_brick,rcl=c(-1,-1,NA),right=NA,
                             filename="mim_survey_2013_multiband_zero.tif",overwrite=T)
names(mim_brick_zero) <- c("gut_stems","til_stems","gut_ext_cm","til_ext_cm")


##Aggregates the data to 9m, 27m and 81m.
mim_brick_9m <- aggregate(mim_brick_zero,fact=3,fun=sum,na.rm=T,filename="mim_survey_2013_9m.tif",
                          overwrite=T)
mim_brick_27m <- aggregate(mim_brick_zero,fact=9,fun=sum,na.rm=T,filename="mim_survey_2013_27m.tif",
                           overwrite=T)
mim_brick_81m <- aggregate(mim_brick_zero,fact=27,fun=sum,na.rm=T,filename="mim_survey_2013_81m.tif")

##Creates a spatial points data frame for all cells that are not NA.
mim_pts_3m <- rasterToPoints(mim_brick_zero,spatial=T,progress="text")
mim_pts_9m <- rasterToPoints(mim_brick_9m,spatial=T,progress="text")
mim_pts_27m <- rasterToPoints(mim_brick_27m,spatial=T,progress="text")
mim_pts_81m <- rasterToPoints(mim_brick_81m,spatial=T,progress="text")

colnames(mim_pts_3m@data) <- c("gut_stems","til_stems","gut_ext_cm","til_ext_cm")
colnames(mim_pts_9m@data) <- c("gut_stems","til_stems","gut_ext_cm","til_ext_cm")
colnames(mim_pts_27m@data) <- c("gut_stems","til_stems","gut_ext_cm","til_ext_cm")
colnames(mim_pts_81m@data) <- c("gut_stems","til_stems","gut_ext_cm","til_ext_cm")

##Creates new columns for presence/absence.
mim_pts_3m@data$til_pres <- mim_pts_3m@data$til_stems
mim_pts_3m@data$til_pres[mim_pts_3m@data$til_pres>0] <- 1
mim_pts_3m@data$gut_pres <- mim_pts_3m@data$gut_stems
mim_pts_3m@data$gut_pres[mim_pts_3m@data$gut_pres>0] <- 1

mim_pts_9m@data$til_pres <- mim_pts_9m@data$til_stems
mim_pts_9m@data$til_pres[mim_pts_9m@data$til_pres>0] <- 1
mim_pts_9m@data$gut_pres <- mim_pts_9m@data$gut_stems
mim_pts_9m@data$gut_pres[mim_pts_9m@data$gut_pres>0] <- 1

mim_pts_27m@data$til_pres <- mim_pts_27m@data$til_stems
mim_pts_27m@data$til_pres[mim_pts_27m@data$til_pres>0] <- 1
mim_pts_27m@data$gut_pres <- mim_pts_27m@data$gut_stems
mim_pts_27m@data$gut_pres[mim_pts_27m@data$gut_pres>0] <- 1

mim_pts_81m@data$til_pres <- mim_pts_81m@data$til_stems
mim_pts_81m@data$til_pres[mim_pts_81m@data$til_pres>0] <- 1
mim_pts_81m@data$gut_pres <- mim_pts_81m@data$gut_stems
mim_pts_81m@data$gut_pres[mim_pts_81m@data$gut_pres>0] <- 1

####Writes the points to shapefile format.
setwd("~/Dropbox/Research/monkeyflower/data/survey_data/rasterized")
writeSpatialShape(mim_pts_3m, "survey_data_gridded_3m", factor2char = TRUE, max_nchar=254)
writeSpatialShape(mim_pts_9m, "survey_data_gridded_9m", factor2char = TRUE, max_nchar=254)
writeSpatialShape(mim_pts_27m, "survey_data_gridded_27m", factor2char = TRUE, max_nchar=254)
writeSpatialShape(mim_pts_81m, "survey_data_gridded_81m", factor2char = TRUE, max_nchar=254)
