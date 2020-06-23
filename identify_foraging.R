# Identify foraging locations
library(ggplot2)
#library(momentuHMM)
library(EMbC)
library(dplyr)
library(sf)

# Pull in trip quality table
t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')

# REMOVE 3 unwanted datasets
t_qual<-t_qual[-which(t_qual$ID%in%c('MEND2_RFBO_Christmas', 'MEND3_RFBO_Christmas',
                                     'AUST3_RFBO_Little Cayman')),]

# make dset summary table

dat_sumr<-t_qual%>%filter(manual_keep=='Y')%>%group_by(ID)%>%
                        summarise(mean_dur=mean(duration, na.rm=T),
                        p_2dur=length(which(duration<2))/n(), n_20loc=length(which(n_locs<21)),
                        mean_cold=mean(max_dist, na.rm=T),
                        p_10col=length(which(max_dist<10))/n())
#write.csv(dat_sumr, 'C:/seabirds/data/dataID_mintrip_decisions.csv', quote=F, row.names=F)
# Booby trial
#t_qual<-t_qual[grep('BRBO', t_qual$ID),]
#qplot(data=t_qual, x=duration, geom='histogram')+facet_wrap(~ID, scales='free')
# BRBO intervals between 1 and 5 min
#t_qual[which(t_qual$duration>72),]$auto_keep<-'N'
#hist(t_qual[which(is.na(t_qual$duration)),]$max_dist) #NAs look ok


for(i in unique(t_qual$ID))
  
{
  cleand<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/clean/',i, '.csv'))
  
 
  # filter out bad trips
  cleand<-cleand[cleand$trip_id%in% t_qual[t_qual$ID==i & t_qual$manual_keep=='Y',]$trip_id,]
  
  # rm troublesome trip
  if(i == 'GILM4_MABO_Clarion') {cleand<-cleand[cleand$trip_id!='MABO_CLR_SDC05_inc_Female_1-11',]}
  
  # classic EMbC 
  cleand$DateTime <- as.POSIXct(strptime(as.character(cleand$DateTime),format="%Y-%m-%d %H:%M:%S",tz="GMT")) # time has to be POSIXct
  forembc <- cleand[,c('DateTime', 'Longitude', 'Latitude', 'trip_id')] #time, longitude, latitude, ID ( ! order is important !)
 
  if(i =='AUST4_MAFR_Little Cayman')
     {BC <- stbc(forembc, spdLim = 1000)}else{
      BC <- stbc(forembc)}
  
  smoothedBC <- smth(BC, dlta=1)
  #view(smoothedBC)

  cleand$embc <- smoothedBC@A
  cleand$embc <- gsub("2","foraging",cleand$embc)
  cleand$embc <- gsub("1","resting",cleand$embc)
  cleand$embc <- gsub("3","commuting",cleand$embc)
  cleand$embc <- gsub("4","relocating",cleand$embc)
  cleand$embc <- gsub("5","DD",cleand$embc)
  
  #write out raw
  write.csv(cleand, paste0('C:/seabirds/sourced_data/tracking_data/foraging_embc/',
                         i, '.csv'), quote=F, row.names=F)
  print(i)
  print(Sys.time())
}

# import csv with embc foraging, export as shapefile for easy QGIS vis

l1<-list.files('C:/seabirds/sourced_data/tracking_data/foraging_embc')
for(i in l1)
{
  p1<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/foraging_embc/',
         i))
 
  p1$ID=i
  p1$ID<-substr(p1$ID, 1, (nchar(p1$ID)-4))
  p1$sp<-do.call(c, lapply(strsplit(as.character(p1$ID), '_'), function(x)x[2]))
  p1<-dplyr::select(p1, ID, sp, trip_id, Latitude, Longitude, DateTime, ColDist, embc)
  p1$trip_id<-as.character(p1$trip_id)
  p1$embc<-as.character(p1$embc)
  p1$DateTime<-as.character(p1$DateTime)
  p1<-p1%>%st_as_sf(coords=c('Longitude', 'Latitude'), crs=4326)
  if(which(i==l1)==1){p2<-p1}else{p2<-rbind(p2, p1)}
  print(i)
}


#st_write(p2, 'C:/seabirds/data/GIS/MASTER_for.shp')

st_write(filter(p2, sp=='BRBO')%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/BRBO_for.shp')
st_write(filter(p2, sp=='MABO')%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/MABO_for.shp')
st_write(filter(p2, sp=='RFBO')%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/RFBO_for.shp')
st_write(filter(p2, sp=='WTSH')%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/WTSH_for.shp')
st_write(filter(p2, sp %in% c('GRFR', 'LEFR', 'MAFR'))%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/FRBD_for.shp')
st_write(filter(p2, sp %in% c('RBTB', 'RTTB'))%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/TRBD_for.shp')
st_write(filter(p2, sp=='SOTE')%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/SOTE_for.shp')
st_write(filter(p2, sp%in% c('BRNO', 'LENO', 'BLNO'))%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/NODD_for.shp')
st_write(filter(p2, sp%in% c('CRTE', 'ROTE', 'CATE'))%>%select(ID, trip_id, embc), 'C:/seabirds/data/GIS/TERN_for.shp')

# not sure if this works **
p2$Longitude<-st_coordinates(p2)[,1]
p2$Latitude<-st_coordinates(p2)[,2]
# **
st_geometry(p2)<-NULL
write.csv(p2, 'C:/seabirds/sourced_data/tracking_data/tracking_master_forage.csv', quote=F, row.names=F)


### OLD ####


# classic EMbC 
cleand$DateTime <- as.POSIXct(strptime(as.character(cleand$DateTime),format="%Y-%m-%d %H:%M:%S",tz="GMT")) # time has to be POSIXct
forembc <- cleand[,c('DateTime', 'Longitude', 'Latitude', 'trip_id')] #time, longitude, latitude, ID ( ! order is important !)
BC <- stbc(forembc)
smoothedBC <- smth(BC, dlta=1)
#view(smoothedBC)

cleand$embc <- smoothedBC@A
cleand$embc <- gsub("2","foraging",cleand$embc)
cleand$embc <- gsub("1","resting",cleand$embc)
cleand$embc <- gsub("3","commuting",cleand$embc)
cleand$embc <- gsub("4","relocating",cleand$embc)
cleand$embc <- gsub("5","DD",cleand$embc)

qplot(data=cleand, x=Longitude, y=Latitude, colour=embc)+facet_wrap(~ID, scales='free')


library(momentuHMM)
### Load raw 
datarawHaggis<-read.csv("rawHaggises.csv")
### Process data
processedHaggis<-prepData(data=rawHaggis,covNames=c("slope","temp"))
### Fit HMM# initial step distribution natural scale 
parametersstepPar0 <- c(1,5,0.5,3)# (mu_1,mu_2,sd_1,sd_2)
# initial angle distribution natural scale 
parametersanglePar0 <- c(0,0,1,8)# (mean_1,mean_2,concentration_1,concentration_2)
fitHaggis <- fitHMM(data = processedHaggis, nbStates = 2,
                    dist = list(step = "gamma", angle = "vm"),
                    Par0 = list(step = stepPar0, angle = anglePar0),
                    formula = ~ slope + I(slope^2),
                    estAngleMean = list(angle=TRUE))



tempdata<-prepData(data.frame(x=d1$longitude, y=d1$latitude,
                              ID=d1$trackID), type="LL")
# could use Colony and ColDist as covariates
plot(tempdata) # this gives idea of starting parameter values
plot(tempdata, compact=T) #better plot

stateNames <- c("resting","foraging", "transiting")

mu0 <- c(0.01,1,2)
sigma0 <- c(0.01,1,2)
zeromass0 <- c(0.5,0.01,0.01)
stepPar0 <- c(mu0,sigma0, zeromass0)
angleMean0 <- c(0,pi,0)
kappa0 <- c(2,1,2)
anglePar0 <- c(angleMean0,kappa0)

fitbird <- fitHMM(data = tempdata, nbStates = 3,
                    dist = list(step = "gamma", angle = "vm"),
                    Par0 = list(step = stepPar0, angle = anglePar0),
                  estAngleMean=list(angle=TRUE),
                   stateNames = stateNames)

mu0 <- c(0.5,2,5)
sigma0 <- c(0.5,2,5)
zeromass0 <- c(0.5,0.1,0.1)
stepPar0 <- c(mu0,sigma0, zeromass0)
angleMean0 <- c(0,pi,0)
kappa0 <- c(2,1,2)
anglePar0 <- c(angleMean0,kappa0)

hmm_1<- fitHMM(data=tempdata,nbStates=3,stepPar0=stepPar0,
               anglePar0=anglePar0)

states<-viterbi(hmm_1)
dat$HMMstates<-states
write.csv(dat, "GPS_141516_clean_resamp_tripsplit_hmm.csv", quote=F, row.names=F)

tempdata<-momentuHMM::prepData(data.frame(x=cleand$Longitude, y=cleand$Latitude,
                                          ID=cleand$trip_id), type="LL")
# could use Colony and ColDist as covariates
#plot(tempdata) # this gives idea of starting parameter values
#plot(tempdata, compact=T) #better plot

print(Sys.time())
stateNames <- c("resting","foraging", "transiting")

mu0 <- c(0.01,0.3,0.6)
sigma0 <- c(0.01,0.3,0.6)
zeromass0 <- c(0.5,0.01,0.01)
stepPar0 <- c(mu0,sigma0, zeromass0)
angleMean0 <- c(0,pi,0)
kappa0 <- c(2,1,2)
anglePar0 <- c(angleMean0,kappa0)

fitbird <- momentuHMM::fitHMM(data = tempdata, nbStates = 3,
                              dist = list(step = "gamma", angle = "vm"),
                              Par0 = list(step = stepPar0, angle = anglePar0),
                              estAngleMean=list(angle=TRUE),
                              stateNames = stateNames)

cleand$HMM<-momentuHMM::viterbi(fitbird)