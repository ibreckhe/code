####Fits a spatial Bayesian glm to the data####

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

##Subsamples the 81m data to get the sample size under control.
samp_size <- round(dim(mim_pts_env_81m@data)[1]*(2/3),digits=0)
set.seed(22)
sampled <- sample(1:dim(mim_pts_env_81m@data)[1],size=samp_size)
mim_pts_sample <- mim_pts_env_81m@data[sampled,]
mim_pts_sample <- mim_pts_sample[complete.cases(mim_pts_sample),]

##Holds the rest of the data for testing.
mim_pts_nosamp <- mim_pts_env_81m@data[which(!sampled%in%rownames(mim_pts_sample)),]

##Gets the coords.
gut_coords <- as.matrix(mim_pts_sample[,c("x","y")])/1000

##Fits a simple glm
y <- mim_pts_sample$gut_pres
log_elev <- scale(log(mim_pts_sample$elev))
log_str_dist <- scale(log(mim_pts_sample$stream_dis+1))
log_can_ht <- scale(log(mim_pts_sample$can_ht+1))

gut_gam_int <- glm(y~log_str_dist+log_elev*log_can_ht,family="binomial")
summary(gut_gam_int)

##spGLM setup parameters.
beta_start <- coef(gut_gam_int)
beta_tune <- t(chol(vcov(gut_gam_int)))
n_batch <- 50
batch_length <- 200
n_samples <- batch_length*n_batch

##Sets initial values for the three chains.
inits <- list()
inits[[1]] <- list(beta = beta_start, phi = runif(1,0.05,20), sigma.sq = runif(1,0.001,0.2), w = 0)
inits[[2]] <- list(beta = beta_start, phi = runif(1,0.05,20), sigma.sq = runif(1,0.001,0.2), w = 0)
inits[[3]] <- list(beta = beta_start, phi = runif(1,0.05,20), sigma.sq = runif(1,0.001,0.2), w = 0)


##does MCMC to estimate model parameters (3 chains in parallel)
##Sets up timing.
start_t <- Sys.time()

##Sets up the MCMC in parallel.
gut_spb_int <- foreach(i=1:3) %dopar% 
  (spGLM(y~log_str_dist+log_elev*log_can_ht,family="binomial",coords=gut_coords,
         starting = inits[[i]], 
         tuning = list(beta = beta_tune, phi = 0.5, sigma.sq = 0.5, w = 0.5), 
         priors = list(beta.Normal = list(rep(0, 5), rep(100, 5)), phi.Unif = c(0.081,20),sigma.sq.IG=c(3,0.5)),
         amcmc = list(n.batch=n_batch,batch.length=batch_length,accept.rate=0.43),
         cov.model="exponential",
         verbose=T,
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
burn_in <- 0.1 * n_samples
sub_samps <- burn_in:n_samples
print(summary(window(gut_spb_int[[1]]$p.beta.theta.samples, start = burn_in)))
print(gelman.diag(samps))
print(spDiag(gut_spb_int[[1]],start=burn_in,thin=5))

##Writes the model object to disk.
setwd("~/GIS/tmp/")
save(gut_spb_int,file="gut_81m_spGLM.Rdata",compress=TRUE)
