####Script to fit the best logistic regression to the presence-absence data and compute maps of predictions####
####Author: Ian Breckheimer
####Date: 12 November 2013.

####Sets up the workspace####
library(raster)
library(maptools)
library(mgcv)

####Read in data####
##Reads in the environmental data.
setwd("/Users/ian/GIS/")

env_names <-c("rough","can_ht","can_pct","elev","isolt","slope","srad","tci_max","tci","snow_melt",
              "ndvi","lidar_int","pond_stream","stream_dist")
env3 <- brick("env_pred_3m.img",crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env9 <- brick("env_pred_9m.img",crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env27 <- brick("env_pred_27m.img",crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
env81 <- brick("env_pred_81m.img",crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

names(env3) <- env_names
names(env9) <- env_names
names(env27) <- env_names
names(env81) <- env_names

##Reads in the survey data.
setwd("~/Dropbox/Research/monkeyflower/data/survey_data/rasterized")
mim_pts_clean_3m <- readShapePoints("survey_data_gridded_3m",
                                    proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_clean_9m <- readShapePoints("survey_data_gridded_9m",
                                    proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_clean_27m <- readShapePoints("survey_data_gridded_27m",
                                     proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_clean_81m <- readShapePoints("survey_data_gridded_81m",
                                     proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))


####Samples the raster covariates at the locations with survey data####
mim_pts_env_3m <- extract(env3,mim_pts_clean_3m,method="simple",sp=T)
mim_pts_env_9m <- extract(env9,mim_pts_clean_9m,method="simple",sp=T)
mim_pts_env_27m <- extract(env27,mim_pts_clean_27m,method="simple",sp=T)
mim_pts_env_81m <- extract(env81,mim_pts_clean_81m,method="simple",sp=T)

##Gets rid of the serious outliers
mim_pts_env_3m <- subset(mim_pts_env_3m,subset=gut_ext_cm < 10000 & til_ext_cm<10000)
mim_pts_env_9m <- subset(mim_pts_env_9m,subset=gut_ext_cm < 10000 & til_ext_cm<10000)
mim_pts_env_27m <- subset(mim_pts_env_27m,subset=gut_ext_cm < 20000 & til_ext_cm<20000)
mim_pts_env_81m <- subset(mim_pts_env_81m,subset=gut_ext_cm < 10000 & til_ext_cm<10000)

####Writes the points to shapefile format.
setwd("~/Dropbox/Research/monkeyflower/data/survey_data/rasterized")
writeSpatialShape(mim_pts_env_3m, "survey_data_gridded_env_3m", factor2char = TRUE, max_nchar=254)
writeSpatialShape(mim_pts_env_9m, "survey_data_gridded_env_9m", factor2char = TRUE, max_nchar=254)
writeSpatialShape(mim_pts_env_27m, "survey_data_gridded_env_27m", factor2char = TRUE, max_nchar=254)
writeSpatialShape(mim_pts_env_81m, "survey_data_gridded_env_81m", factor2char = TRUE, max_nchar=254)


####Builds the models####
gut_gam_3m <- glm(gut_pres~I(log(stream_dist+0.01))+elev+can_ht+I(elev^2)+I(rough^2),
                  data=mim_pts_env_3m@data,family="binomial")
gut_gam_9m <- glm(gut_pres~I(log(stream_dist+0.01))+elev+can_ht+I(elev^2)+I(rough^2),
                  data=mim_pts_env_9m@data,family="binomial")
gut_gam_27m <- glm(gut_pres~I(log(stream_dist+0.01))+elev+can_ht+I(elev^2)+I(rough^2),
                   data=mim_pts_env_27m@data,family="binomial")
gut_gam_81m <- glm(gut_pres~I(log(stream_dist+0.01))+elev+can_ht+I(elev^2)+I(rough^2),
                   data=mim_pts_env_81m@data,family="binomial")

til_gam1_3m <- glm(til_pres~I(log(stream_dist+0.01))+elev+can_ht+rough+I(rough^2)+I(can_ht^2)+I(log(stream_dist+0.01)):I(elev^2),
                   data=mim_pts_env_3m@data,family="binomial")
til_gam1_9m <- glm(til_pres~I(log(stream_dist+0.01))+elev+can_ht+rough+I(rough^2)+I(can_ht^2)+I(log(stream_dist+0.01)):I(elev^2),
                   data=mim_pts_env_9m@data,family="binomial")
til_gam1_27m <- glm(til_pres~I(log(stream_dist+0.01))+elev+can_ht+rough+I(rough^2)+I(can_ht^2)+I(log(stream_dist+0.01)):I(elev^2),
                    data=mim_pts_env_27m@data,family="binomial")
til_gam1_81m <- glm(til_pres~I(log(stream_dist+0.01))+elev+can_ht+rough+I(rough^2)+I(can_ht^2)+I(log(stream_dist+0.01)):I(elev^2),
                    data=mim_pts_env_81m@data,family="binomial")

model_list <- list(gut_gam_3m,gut_gam_9m,gut_gam_27m,gut_gam_81m,til_gam1_3m,til_gam1_9m,til_gam1_27m,til_gam1_81m)


###Assesses each model using AUC####
model_list <- list(gut_gam_3m,gut_gam_9m,gut_gam_27m,gut_gam_81m,til_gam1_3m,til_gam1_9m,til_gam1_27m,til_gam1_81m)
mod_names <- c("gut_gam_3m","gut_gam_9m","gut_gam_27m","gut_gam_81m","til_gam1_3m","til_gam1_9m","til_gam1_27m","til_gam1_81m")
data_list <- list(mim_pts_env_3m@data,mim_pts_env_9m@data,mim_pts_env_27m@data,mim_pts_env_81m@data,
                  mim_pts_env_3m@data,mim_pts_env_9m@data,mim_pts_env_27m@data,mim_pts_env_81m@data)

gut_evals <- as.list(rep(0,4))
for (i in seq(1:4)){
  flush.console()
  print(paste("Evaluating ",mod_names[i],sep=""))
  
  presence <- data_list[[i]][data_list[[i]]$gut_pres == 1,]
  absence <- data_list[[i]][data_list[[i]]$gut_pres == 0,]
  
  ##Creates a ROC plot.
  gut_evals[[i]] <- evaluate(presence, absence, model_list[[i]])
}

til_evals <- as.list(rep(0,4))
for (i in seq(5:8)){
  flush.console()
  print(paste("Evaluating ",mod_names[i],sep=""))
  presence <- data_list[[i]][data_list[[i]]$til_pres == 1,]
  absence <- data_list[[i]][data_list[[i]]$til_pres == 0,]
  
  ##Creates a ROC plot.
  til_evals[[i]] <- evaluate(presence, absence, model_list[[i]])
}

##Plots the ROC plots.
par(mfrow=c(4,2),mar=c(1,1,0,0),oma=c(4,4,2,1))
for (i in 1:4){
  ##tilingii plots
  plot(til_evals[[i]],'ROC',ann=F,axes=F)
  box()
  text(0.6,0.2,labels=paste("AUC=",round(til_evals[[i]]@auc,3),sep=""))
  axis(2,labels=T,tick=T)
  if(i == 4) axis(1,tick=T) 
  else axis(side=1,labels=F,tick=T)
  if(i == 1) title(main="tilingii",font=3,outer=F,xpd=NA,line=1) 
  ##guttatus plots
  plot(gut_evals[[i]],'ROC',ann=F,axes=F)
  box()
  text(0.6,0.2,labels=paste("AUC=",round(gut_evals[[i]]@auc,3),sep=""))
  axis(2,labels=F,tick=T)
  if(i == 4) axis(1,labels=T,tick=T) 
  else axis(side=1,labels=F,tick=T)
  if(i == 1) title(main="guttatus",font=3,outer=F,xpd=NA,line=1) 
}
mtext(text="False Positive Rate",side=1,outer=T,padj=2)
mtext(text="True Positive Rate",side=2,outer=T,padj=-2)


####Returns maps of predictions and standard errors for each model####
predfun <- function(model, data) {
  v <- predict(model, data, se.fit=TRUE)
  cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
}

gam1_pred3 <- predict(env3,model=gut_gam_3m,overwrite=T,fun=predfun,index=1:2,progress="text")
gam1_pred9 <- predict(env9,model=gut_gam_9m,overwrite=T,fun=predfun,index=1:2,progress="text")
gam1_pred27 <- predict(env27,model=gut_gam_27m,overwrite=T,fun=predfun,index=1:2,progress="text")
gam1_pred81 <- predict(env81,model=gut_gam_81m,overwrite=T,fun=predfun,index=1:2,progress="text")

til_gam1_pred3 <- predict(env3,model=til_gam1_3m,overwrite=T,fun=predfun,index=1:2,progress="text")
til_gam1_pred9 <- predict(env9,model=til_gam1_9m,overwrite=T,fun=predfun,index=1:2,progress="text")
til_gam1_pred27 <- predict(env27,model=til_gam1_27m,overwrite=T,fun=predfun,index=1:2,progress="text")
til_gam1_pred81 <- predict(env81,model=til_gam1_81m,overwrite=T,fun=predfun,index=1:2,progress="text")

pred_list <- list(gam1_pred81,gam1_pred27,gam1_pred9,gam1_pred3,
                  til_gam1_pred81,til_gam1_pred27,til_gam1_pred9,til_gam1_pred3)
pred_names <- c("gam1_pred81","gam1_pred27","gam1_pred9","gam1_pred3",
                "til_gam1_pred81","til_gam1_pred27","til_gam1_pred9","til_gam1_pred3")
for (i in 1:length(pred_list)) {
  names(pred_list[[i]])[1] <- pred_names[i]
}

##Adds 95% confidence intervals to the stacks of predictions.
setwd("/Users/ian/GIS/")
critval <- 1.96 ## approx 95% CI
inv_logit <- function(x) {exp(x[1])/(1+exp(x[1]))}
upr_fun <- function(x){inv_logit(x[1] + critval*x[2])}
lwr_fun <- function(x){inv_logit(x[1] - critval*x[2])}

for (i in 1:length(pred_list)) {
  flush.console()
  print(paste("Now processing raster ",names(pred_list[[i]])[1],sep=""))
  est <- calc(pred_list[[i]],fun=inv_logit,progress="text",overwrite=T)
  upr <- calc(pred_list[[i]],fun=upr_fun,progress="text",overwrite=T)
  lwr <- calc(pred_list[[i]],fun=lwr_fun,progress="text",overwrite=T)
  pred <- brick(stack(est,lwr,upr),filename=paste(names(pred_list[[i]])[1],".tif",sep=""),overwrite=T,bylayer=F)
}

####Assesses how prevalence changes across scales####

##Total area surveyed.
tot_area <- (length(gut_evals[[4]]@presence)+length(gut_evals[[4]]@absence))*81^2
scales <- c(3,9,27,81)
cell_area <- scales^2

gut_prev <- c()
til_prev <- c()
for (i in 1:length(scales)){
  gut_prev <- c(gut_prev,(length(gut_evals[[i]]@presence)*cell_area[i])/tot_area)
  til_prev <- c(til_prev,(length(til_evals[[i]]@presence)*cell_area[i])/tot_area)
}

prev_scale <- data.frame(scales,gut_prev,til_prev)

##Fits a power-law model.
gut_scale_lm <- lm(I(log(gut_prev))~I(log(scales)),data=prev_scale)
til_scale_lm <- lm(I(log(til_prev))~I(log(scales)),data=prev_scale)

##Gets model predictions for 5cm scale.
gut_scale_area <- predict(gut_scale_lm,newdata=data.frame(scales=c(0.0005,0.1,1)),se.fit=T,interval="prediction",
                          type="response")
til_scale_area <- predict(til_scale_lm,newdata=data.frame(scales=c(0.0005,0.1,1)),se.fit=T,interval="prediction",
                          type="response")
gut_est_area_m <- tot_area*exp(gut_scale_area$fit[1,1])
til_est_area_m <- tot_area*exp(til_scale_area$fit[1,1])

##Compares to measured area.
gut_meas_area_m <- sum(mim_pts_env_3m$gut_ext_cm)*0.0001
til_meas_area_m <- sum(mim_pts_env_3m$til_ext_cm)*0.0001
gut_meas_prev <- gut_meas_area_m/tot_area
til_meas_prev <- til_meas_area_m/tot_area

####Plots outputs.####
setwd("/Users/ian/Dropbox/Research/monkeyflower/plots")

##Plots the models
svg("grain_prevalence.svg",width=5,height=5)

gut_col <- rgb(0.4,0.9,0.4,0.5)
til_col <- rgb(0.4,0.4,0.9,0.5)

plot(log(gut_prev)~log(scales),data=prev_scale,col=gut_col,pch=20,xlim=c(-6,5),ylim=c(-20,0),
     xlab="log(Grain Size (m))",ylab="log(Prevalence)",axes=F,cex=1.5)
axis(1,labels=as.character(c(0.0001,0.001,0.01,0.1,1,10,100)),at=log(c(0.0001,0.001,0.01,0.1,1,10,100)))
axis(2,labels=as.character(c(1e-9,1e-7,1e-5,1e-3,1e-1,1)),at=log(c(1e-9,1e-7,1e-5,1e-3,1e-1,1)))
points(log(til_prev)~log(scales),data=prev_scale,col=til_col,pch=20,cex=1.5)
abline(gut_scale_lm,lty=2,col=gut_col)
abline(til_scale_lm,lty=2,col=til_col)
legend("top",legend=c("guttatus","tilingii"),bty="n",pch=c(20,20),col=c(gut_col,til_col),horiz=T,cex=1.5)

##Puts the field estimates on a plot.
points(log(0.0025),log(gut_meas_prev),pch=21,col=gut_col,cex=1.5)
points(log(0.0025),log(til_meas_prev),pch=21,col=til_col,cex=1.5)

dev.off()
