####Script to analyze spatial autocorrelation in the survey data.
####Author: Ian Breckheimer
####Date:10 December 2013

####Sets up the workspace####
library(geoR)
library(maptools)
setwd("~/Dropbox/Research/monkeyflower/data/survey_data/rasterized/")

####Brings in the data.
mim_pts_clean_3m <- readShapePoints("survey_data_gridded_env_3m",
                                    proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_clean_9m <- readShapePoints("survey_data_gridded_env_9m",
                                    proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_clean_27m <- readShapePoints("survey_data_gridded_env_27m",
                                     proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_clean_81m <- readShapePoints("survey_data_gridded_env_81m",
                                     proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))

####Analyzes spatial autocorrelation in the abundance data
##Pulls out just the cells with presence data.
gut_stems_nonzero <- mim_pts_clean_81m[mim_pts_clean_81m$gut_pres==1,]
til_stems_nonzero <- mim_pts_clean_81m[mim_pts_clean_81m$til_pres==1,]

##Computes the empirical variogram for plots with at least one stem.
gut_stems_vario <- variog(coords=gut_stems_nonzero@coords,
                          data=log(gut_stems_nonzero@data$gut_stems),
                          breaks=seq(0,4000,by=30),
                          option="bin")
til_stems_vario <- variog(coords=til_stems_nonzero@coords,
                          data=log(til_stems_nonzero@data$til_stems),
                          breaks=seq(0,4000,by=30),
                          option="bin")

##Computes the variogram for the presence/absence data.
gut_pres_vario <- variog(coords=mim_pts_clean_81m@coords,
                          data=mim_pts_clean_81m@data$gut_pres,
                          breaks=seq(0,4000,by=30),
                          option="bin")
til_pres_vario <- variog(coords=mim_pts_clean_81m@coords,
                          data=mim_pts_clean_81m@data$til_pres,
                          breaks=seq(0,4000,by=30),
                          option="bin")


####Writes plots####
setwd("~/Dropbox/Research/monkeyflower/plots/")

##Colors for the two species
gut_col <- rgb(0.4,0.9,0.4,0.5)
til_col <- rgb(0.4,0.4,0.9,0.5)

##Plots the counts and log counts.
svg("abundance_hist.svg",width=5,height=5)
par(mfrow=c(2,1),oma=c(4,4,2,2),mar=c(1,1,1,1),xpd=F)
hist(log(gut_stems_nonzero$gut_stems),breaks=seq(0,12,by=1),main="",axes=F,xlab="",col=gut_col)
axis(1,labels=F,at=seq(0,12,by=2))
axis(2,labels=T)
hist(log(til_stems_nonzero$til_stems),breaks=seq(0,12,by=1),main="",col=til_col)
mtext(1,text="Log Abundance",outer=T,line=2)
mtext(2,text="Frequency",outer=T,line=2)
legend("topright",legend=c("guttatus","tilingii"),bty="n",fill=c(gut_col,til_col))
dev.off()


##Plots the variogram for abundance.
svg("abundance_variog.svg",width=5,height=5)
par(mfrow=c(2,1),oma=c(4,4,2,2),mar=c(1,1,1,1),xpd=F)
plot(gut_stems_vario,pch=20,col=gut_col,axes=F,ylim=c(0,10))
axis(1,labels=F,at=seq(0,4000,by=1000))
axis(2,labels=T)

plot(til_stems_vario,pch=20,col=til_col,axes=F,ylim=c(0,10))
axis(1,labels=T,ticks=T)
axis(2,labels=T,ticks=T)
mtext(1,text="lag distance (m)",outer=T,line=2)
mtext(2,text="semivariance",outer=T,line=2)
legend("topright",legend=c("guttatus","tilingii"),bty="n",pch=c(20,20),col=c(gut_col,til_col),horiz=T)
dev.off()

##Plots the variogram for abundance.
svg("presence_variog.svg",width=5,height=5)
par(mfrow=c(2,1),oma=c(4,4,2,2),mar=c(1,1,1,1),xpd=F)
plot(gut_pres_vario,pch=20,col=gut_col,axes=F,ylim=c(0,0.2))
axis(1,labels=F,at=seq(0,4000,by=1000))
axis(2,labels=T)

plot(til_pres_vario,pch=20,col=til_col,axes=F,ylim=c(0,0.2))
axis(1,labels=T,ticks=T)
axis(2,labels=T,ticks=T)
mtext(1,text="lag distance (m)",outer=T,line=2)
mtext(2,text="semivariance",outer=T,line=2)
legend("topright",legend=c("guttatus","tilingii"),bty="n",pch=c(20,20),col=c(gut_col,til_col),horiz=T)
dev.off()