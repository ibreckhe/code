####Fits a spatial Bayesian glm to the gridded survey data####
####Author:Ian Breckheimer.
####Date:13 December 2013

####Sets up the workspace####
library(spBayes)
library(maptools)
library(doParallel)
registerDoParallel(cores=3)

##Reads in the survey data.
setwd("~/Dropbox/Research/monkeyflower/data/survey_data/rasterized")
mim_pts_env_3m <- readShapePoints("survey_data_gridded_env_3m",
                                  proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_env_9m <- readShapePoints("survey_data_gridded_env_9m",
                                  proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_env_27m <- readShapePoints("survey_data_gridded_env_27m",
                                   proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
mim_pts_env_81m <- readShapePoints("survey_data_gridded_env_81m",
                                   proj4string=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))

##Subsamples the 81m data to get the sample size under control.(1/8 of the data=3min,2/3 of the data=4hr)
samp_size <- round(dim(mim_pts_env_81m@data)[1]*(8/8),digits=0)
set.seed(22)
sampled <- sample(1:dim(mim_pts_env_81m@data)[1],size=samp_size)
mim_pts_sample <- mim_pts_env_81m@data[sampled,]
mim_pts_sample <- mim_pts_sample[complete.cases(mim_pts_sample),]

##Holds the rest of the data for testing.
mim_pts_nosamp <- mim_pts_env_81m@data[which(!sampled%in%rownames(mim_pts_sample)),]

##Gets the coords and converts to kilometers to avoid numerical problems
obs_coords <- as.matrix(mim_pts_sample[,c("x","y")])/1000

##Centers and scales predictors
gut_pres <- mim_pts_sample$gut_pres
til_pres <- mim_pts_sample$til_pres
log_elev <- scale(log(mim_pts_sample$elev))
log_str_dist <- scale(log(mim_pts_sample$stream_dis+1))
log_can_ht <- scale(log(mim_pts_sample$can_ht+1))
srad <- scale(mim_pts_sample$srad)

##Fits the best simple glm for both species.
gut_glm_best <- glm(gut_pres~log_str_dist+log_elev*log_can_ht,family="binomial")
summary(gut_glm_best)

til_glm_best <- glm(til_pres~log_str_dist+srad,
                    family="binomial")
summary(til_glm_best)

####Fits the guttatus spatial model####

##spGLM setup parameters.
beta_start <- coef(gut_glm_best)
beta_tune <- t(chol(vcov(gut_glm_best)))
n_batch <- 50
batch_length <- 200
n_samples <- batch_length*n_batch

##Sets initial values for the three chains.
inits <- list()
inits[[1]] <- list(beta = beta_start, phi = runif(1,0.081,10), sigma.sq = runif(1,0.001,0.2), w = 0)
inits[[2]] <- list(beta = beta_start, phi = runif(1,0.081,10), sigma.sq = runif(1,0.001,0.2), w = 0)
inits[[3]] <- list(beta = beta_start, phi = runif(1,0.081,10), sigma.sq = runif(1,0.001,0.2), w = 0)

##does MCMC to estimate model parameters (3 chains in parallel)
##Sets up timing.
start_t <- Sys.time()

##Sets up the MCMC in parallel.
gut_spb_int <- foreach(i=1:3) %dopar% 
  (spGLM(gut_pres~log_str_dist+log_elev*log_can_ht,family="binomial",coords=obs_coords,
         starting = inits[[i]], 
         tuning = list(beta = beta_tune, phi = 0.5, sigma.sq = 0.1, w = 0.5), 
         priors = list(beta.Normal = list(rep(0, 5), rep(100, 5)), phi.Unif = c(0.081,10),sigma.sq.IG=c(6,0.5)),
         cov.model="exponential",
         verbose=T,
         amcmc = list(n.batch=n_batch,batch.length=batch_length,accept.rate=0.43),
         n.report=50))

##Reports the elapsed time.
end_t <- Sys.time()
elapsed <- end_t - start_t
print(paste("Total Execution Time:",round(elapsed,2),units(elapsed)))

##Plots diagnostics of the chains.
samps <- mcmc.list(gut_spb_int[[1]]$p.beta.theta.samples,
                   gut_spb_int[[2]]$p.beta.theta.samples,
                   gut_spb_int[[3]]$p.beta.theta.samples)
plot(samps)

##Checks model outputs
burn_in <- 0.5 * n_samples
sub_samps <- burn_in:n_samples
print(summary(window(gut_spb_int[[1]]$p.beta.theta.samples, start = burn_in)))
print(gelman.diag(samps))
print(spDiag(gut_spb_int[[1]],start=burn_in,thin=5))

##Compares credible intervals with nonspatial frequentist standard errors.
summary(gut_glm_best)

##Writes the model object to disk.
setwd("~/GIS/tmp/")
save(gut_spb_int,file="gut_81m_spGLM.Rdata",compress=TRUE)

####Fits the tilingii model in a spatial Bayesian framework.####

##spGLM setup parameters.
til_beta_start <- coef(til_glm_best)
til_beta_tune <- t(chol(vcov(til_glm_best)))
n_batch <- 50
batch_length <- 200
n_samples <- batch_length*n_batch

##Sets initial values for the three chains.
til_inits <- list()
til_inits[[1]] <- list(beta = til_beta_start, phi = runif(1,0.05,5), sigma.sq = runif(1,0.04,0.06), w = 0)
til_inits[[2]] <- list(beta = til_beta_start, phi = runif(1,0.05,5), sigma.sq = runif(1,0.04,0.06), w = 0)
til_inits[[3]] <- list(beta = til_beta_start, phi = runif(1,0.05,5), sigma.sq = runif(1,0.04,0.06), w = 0)

##Sets up timing.
start_t <- Sys.time()

##Sets up the MCMC in parallel.
til_spb_int <- foreach(i=1:3) %dopar% 
  (spGLM(til_pres~log_str_dist+srad,family="binomial",coords=obs_coords,
         starting = til_inits[[i]], 
         tuning = list(beta = til_beta_tune, phi = 0.5, sigma.sq = 0.1, w = 0.5), 
         priors = list(beta.Normal = list(rep(0, 3), rep(100, 3)), phi.Unif = c(0.081,10),sigma.sq.IG=c(5,0.5)),
         cov.model="exponential",
         verbose=T,
         amcmc = list(n.batch=n_batch,batch.length=batch_length,accept.rate=0.43),
         n.report=50))

##Reports the elapsed time.
end_t <- Sys.time()
elapsed <- end_t - start_t
print(paste("Total Execution Time:",round(elapsed,2),units(elapsed)))

##Plots diagnostics of the chains.
til_samps <- mcmc.list(til_spb_int[[1]]$p.beta.theta.samples,
                   til_spb_int[[2]]$p.beta.theta.samples,
                   til_spb_int[[3]]$p.beta.theta.samples)
plot(til_samps)

##Checks model outputs
burn_in <- 0.5 * n_samples
sub_samps <- burn_in:n_samples
print(summary(window(til_spb_int[[1]]$p.beta.theta.samples, start = burn_in,thin=5)))
print(gelman.diag(til_samps))
print(spDiag(til_spb_int[[1]],start=burn_in,thin=5))

##Compares credible intervals with nonspatial frequentist standard errors.
summary(til_glm_best)

##Writes the model object to disk.
save(til_spb_int,file="til_81m_spGLM.Rdata",compress=TRUE)


####Creates a model for abundance####

##Extracts data with at least one guttatus stem
gut_stems_nonzero <- mim_pts_env_27m[mim_pts_env_27m$gut_pres==1,]
gut_stems_nonzero$log_gut_stems <- log(gut_stems_nonzero$gut_stems)

##Scales covariates
gut_abund <- gut_stems_nonzero$log_gut_stems
log_elev <- scale(log(gut_stems_nonzero$elev))
log_str_dist <- scale(log(gut_stems_nonzero$stream_dis+1))
log_can_ht <- scale(log(gut_stems_nonzero$can_ht+1))
lidar_int <- scale(gut_stems_nonzero$lidar_int)
rough <- scale(gut_stems_nonzero$rough)
srad <- scale(gut_stems_nonzero$srad)
log_slope <- scale(log(gut_stems_nonzero$slope+1))
log_tci_max <- scale(log(gut_stems_nonzero$tci_max))
pond_stream <- scale(gut_stems_nonzero$pond_strea)

##Examines covariates
pairs(data.frame(gut_abund,log_elev,log_str_dist,log_can_ht,rough,srad,log_tci_max,log_slope))

##Builds a reasonable model.
gut_abund_mod <- glm(gut_abund~log_can_ht,family="gaussian")
summary(gut_abund_mod)

##Extracts data with at least one tilingii stem
til_stems_nonzero <- mim_pts_env_27m[mim_pts_env_27m$til_pres==1,]
til_stems_nonzero$log_til_stems <- log(til_stems_nonzero$til_stems)

##Builds a reasonable model.
til_abund_mod <- glm(til_abund~log_elev)
summary(til_abund_mod)
