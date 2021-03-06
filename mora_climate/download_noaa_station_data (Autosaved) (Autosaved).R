##Script to play around with downloading climate data from the world bank.
library(devtools)
library(rnoaa)
library(plyr)
library(reshape2)

##Sets the API key for downloading NOAA data.
options(noaakey="NjuEQnGAERfxLrZNIiHqMUJsSeqiPuRR")

##Sets the working directory.
setwd("~/code/mora_climate")

##Current stations in the vicinity of Mt. Rainier.
# stationIDs <- c("GHCND:USS0021C38S",
#                 "GHCND:USS0021B63S",
#                 "GHCND:USS0021B13S",
#                 "GHCND:USR0000WGRN",
#                 "GHCND:USS0021B62S",
#                 "GHCND:USC00454764",
#                 "GHCND:US1WALW001",
#                 "GHCND:USS0021C17S",
#                 "GHCND:USS0021C40S",
#                 "GHCND:USC00455704",
#                 "GHCND:USC00456262",
#                 "GHCND:USS0021C35S",
#                 "GHCND:USS0021C33S",
#                 "GHCND:USC00456896",
#                 "GHCND:USC00456898",
#                 "GHCND:USC00456909",
#                 "GHCND:USS0021C28S")

##All current stations within 50km of the summit.
stationIDs <- c("GHCND:USS0021C35S",
				"GHCND:USS0021B62S",
				"GHCND:USC00455704",
				"GHCND:USS0021B63S",
				"GHCND:USC00456896",
				"GHCND:USS0021C28S",
				"GHCND:USC00455224",
				"GHCND:USC00456909",
				"GHCND:USS0021C33S",
				"GHCND:USC00456898",
				"GHCND:USS0021C40S",
				"GHCND:USC00456262",
				"GHCND:USS0021B13S",
				"GHCND:USC00454764",
				"GHCND:USS0021C17S",
				"GHCND:USS0021C38S",
				"GHCND:US1WALW0015",
				"GHCND:USS0021B42S",
				"GHCND:USR0000WGRN",
				"GHCND:US1WAPR0028")

# ##All (current and historical) stations in the vicinity of Mt. Rainier.
# stationIDs <- c("GHCND:USC00450285",
                # "GHCND:USC00450945",
                # "GHCND:USC00450969",
                # "GHCND:USS0021C38S",
                # "GHCND:USS0021B63S",
                # "GHCND:USC00451113",
                # "GHCND:USC00451110",
                # "GHCND:USS0021B13S",
                # "GHCND:USC00452493",
                # "GHCND:USC00452722",
                # "GHCND:USC00453177",
                # "GHCND:USC00454286",
                # "GHCND:USC00453219",
                # "GHCND:USR0000WGRN",
                # "GHCND:USC00453357",
                # "GHCND:USR0000WHAG",
                # "GHCND:USS0021B62S",
                # "GHCND:USC00454620",
                # "GHCND:USC00454619",
                # "GHCND:USC00456894",
                # "GHCND:USC00454764",
                # "GHCND:US1WALW0015",
                # "GHCND:USC00455425",
                # "GHCND:USS0021C17S",
                # "GHCND:USS0021C40S",
                # "GHCND:USC00455704",
                # "GHCND:USC00456201",
                # "GHCND:USC00456265",
                # "GHCND:USC00456262",
                # "GHCND:USS0021C35S",
                # "GHCND:USC00456381",
                # "GHCND:USS0021C33S",
                # "GHCND:USC00456892",
                # "GHCND:USC00456896",
                # "GHCND:USC00456898",
                # "GHCND:USC00456900",
                # "GHCND:USC00456909",
                # "GHCND:USC00458653",
                # "GHCND:USS0021C28S",
                # "GHCND:USC00459170",
                # "GHCND:USC00456385",
                # "GHCND:USC00459171")

##Function to get all daily temperature and precipitation data for the stations. 
##NOAA limits API calls to 1000 records, so this is a workaround.
get_ghcn_daily <- function(start_date,end_date,type,stations,chunk_days=20){
  dates <- seq(start_date,end_date,by=chunk_days)
  data_list <- list()
  for (i in 2:length(dates)){
    data_list[[i-1]] <- noaa(datasetid = "GHCND",datatypeid=type,stationid = stations,
                             startdate=as.character(dates[i-1]),enddate=as.character(dates[i]),limit=1000)
  }
  data <- data_list[[1]]
  for (i in 2:length(data_list)){
    data <- noaa_combine(data,data_list[[i]])
  }
  return(data)
}

##Gets the data.
start <- as.Date("2013-10-01")
end <- as.Date("2016-09-30")
tmax_data <- get_ghcn_daily(start,end,type="TMAX",stations=stationIDs,chunk_days=20)[[1]]
tmin_data <- get_ghcn_daily(start,end,type="TMIN",stations=stationIDs,chunk_days=20)[[1]]
prcp_data <- get_ghcn_daily(start,end,type="PRCP",stations=stationIDs,chunk_days=20)[[1]]

##Quick plot to check validity.
par(mfrow=c(3,1))
scatter.smooth(as.Date(tmax_data$date),tmax_data$value,pch=20,cex=0.2,
               col=rgb(0.5,0.5,0.5,0.5),main="TMAX",lpars=list(col="red"),span=0.1)
scatter.smooth(as.Date(tmin_data$date),tmin_data$value,pch=20,cex=0.2,
               col=rgb(0.5,0.5,0.5,0.5),main="TMIN",lpars=list(col="red"),span=0.1)
scatter.smooth(as.Date(prcp_data$date),prcp_data$value,pch=20,cex=0.2,
               col=rgb(0.5,0.5,0.5,0.5),main="PRCP",lpars=list(col="red"),span=0.05)


##Writes the data to disk.
write.csv(tmax_data,"MORA_daily_climate_tmax.csv")
write.csv(tmin_data,"MORA_daily_climate_tmin.csv")
write.csv(prcp_data,"MORA_daily_climate_prcp.csv")

##Brings in the historical station data.
statdat <- read.csv("station_data_50km_1981_2013.csv")

tmax_dat <- statdat[,c("STATION","LATITUDE","LONGITUDE","DATE","TMAX")]
tmax_dat$DATE <- strptime(as.character(tmax_dat$DATE),format="%Y%m%d")

# ##Examines the data.
# scatter.smooth(as.Date(tmax_dat$DATE),tmax_dat$TMAX,pch=20,cex=0.2,
               # col=rgb(0.5,0.5,0.5,0.5),main="TMAX",lpars=list(col="red"))

##Eliminates impossible values
tmax_dat$TMAX[which(tmax_dat$TMAX > 500)] <- NA
tmax_dat$TMAX[which(tmax_dat$TMAX < -300)] <- NA

##Creates a continuous sequence of days.
start <- min(tmax_dat$DATE)
end <- max(tmax_dat$DATE)
days <- data.frame(date=seq(start,end,by="day"))

##Casts the tmax data to "wide" format.
tmax_wide<- dcast(tmax_dat,DATE~STATION,value.var="TMAX")
dim(tmax_wide)

##Computes a median across all stations for each day.
tmax_median <- apply(tmax_wide[,2:dim(tmax_wide)[2]],MARGIN=1,FUN=median,na.rm=T)

##Eliminates all values more than 20C from the median.
for (i in 2:dim(tmax_wide)[2]){
	diff <- abs(tmax_wide[,i] - tmax_median)
	out <- which(diff>150)
	tmax_wide[out,i] <- NA
}

##Collapses data down for plotting.
attributes(tmax_wide$DATE) <- NULL
tmax_nar <- melt(tmax_wide,id.vars="DATE_FAC")
