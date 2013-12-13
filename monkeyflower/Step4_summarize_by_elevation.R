####Script to summarize the GLM habitat model predictions by elevation####
####Author: Ian Breckheimer
####Date: 1 December 2013

####Sets up the workspace
library(raster)
library(doParallel)
registerDoParallel(cores=3)
setwd("~/GIS/")

####Brings in the environmental data####
envs <- c("env_pred_81m.img","env_pred_27m.img","env_pred_9m.img","env_pred_3m.img",
          "env_pred_81m.img","env_pred_27m.img","env_pred_9m.img","env_pred_3m.img")

env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")

env_list <- list()
for(i in 1:length(envs)){
  env_list[[i]] <- brick(envs[i],crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  names(env_list[[i]]) <- env_names
}

####Brings in the predictions and standard errors.
prds <- c("gam1_pred81.img","gam1_pred27.img","gam1_pred9.tif","gam1_pred3.tif"
          "til_gam1_pred81.tif","til_gam1_pred27.tif","til_gam1_pred9.tif","gam1_pred3.tif")

prd_list <- list()
for(i in 1:length(prds)){
  prd_list[[i]] <- brick(prds[i],crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
}

####Plots the relationship between habitat area and elevation.####

##Functions to use in raster calculations.
round_fun <- function(x){round(x/100) * 100}
bound_fun <- function(x){
  bound <- seq(400,4500,by=100)
  return(rep(x,length(bound))==bound)
}

##Function to summarize by elevation.
sum_elev <- function(prediction_brick,elev_raster){
  flush.console()
  print(paste("Processing raster: ", names(prediction_brick)[1], sep=""))
  elev_band <- calc(elev_raster,fun=round_fun,progress="text")
  elev_bands <- calc(elev_band,fun=bound_fun,forceapply=T)
  prd_bands <- prediction_brick[[1]]*elev_bands
  lwr_bands <- prediction_brick[[2]]*elev_bands
  upr_bands <- prediction_brick[[3]]*elev_bands
  res <- (elev_band@extent@xmax - elev_band@extent@xmin) / elev_band@ncols
  prd_sums <- cellStats(prd_bands,stat='sum') * res^2
  lwr_sums <- cellStats(lwr_bands,stat='sum') * res^2
  upr_sums <- cellStats(upr_bands,stat='sum') * res^2
  elev_sums <- cellStats(elev_bands,stat='sum') * res^2
  sum <- data.frame(elev=seq(400,4500,by=100),res=res,model=names(prediction_brick)[1],prd_sums,lwr_sums,upr_sums,available=elev_sums)
  return(sum)
}

##Processes the rasters in parallel.
sums <- foreach(i=1:length(prds),.combine=rbind) %dopar% (sum_elev(prd_list[[i]],env_list[[i]][["elev"]])) 
  
##Converts to ha
sums[,4:7] <- sums[,4:7] / 10000
gut_sums_81m <- sums[sums$model=="gam1_pred81.1",]
til_sums_81m <- sums[sums$model=="til_gam1_pred81.1",]

##Makes the plot.
svg("~/Dropbox/Research/monkeyflower/plots/mim_elev_glm.svg",width=5,height=5)
par(mfrow=c(1,1))
plot(upr_sums~elev,data=til_sums_81m,type="n",ylab="Habitat Area (ha)",xlab="Elevation (m)",xlim=c(400,3000))
polygon(c(gut_sums_81m$elev,rev(gut_sums_81m$elev)),c(gut_sums_81m$lwr_sums,rev(gut_sums_81m$upr_sums)), 
        col=rgb(0.4,0.9,0.4,0.5),border=NULL,lty="blank")
lines(gut_sums_81m$elev,gut_sums_81m$prd_sums,lty=2,col=rgb(0.4,0.9,0.4,1),lwd=2)
polygon(c(til_sums_81m$elev,rev(til_sums_81m$elev)),c(til_sums_81m$lwr_sums,rev(til_sums_81m$upr_sums)), 
        col=rgb(0.4,0.4,0.9,0.5),border=NULL,lty="blank")
lines(til_sums_81m$elev,til_sums_81m$prd_sums,lty=2,col=rgb(0.4,0.4,0.9,1),lwd=2)
legend(2000,700,legend=c("guttatus","tilingii"),bty="n",
       lwd=2,border=0,lty=2,seg.len=1,col=c(rgb(0.4,0.9,0.4,1),rgb(0.4,0.4,0.9,1)),
       fill=c(rgb(0.4,0.9,0.4,0.5),rgb(0.4,0.4,0.9,0.5)))
dev.off()