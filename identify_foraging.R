# Identify foraging locations
library(ggplot2)
#library(momentuHMM)
library(EMbC)
library(dplyr)
library(sf)
# find scale of interaction BL code
source('C:/seabirds/sourced_data/code/new_mIBA/findScaleUPDATEjun20.R')

# Pull in trip quality table
t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')

# Fix yoda dates
#t_qual$departure<-as.character(t_qual$departure)
#t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$departure<-as.character(format(as.Date(as.POSIXlt(as.numeric(
#  t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$departure),
#  origin="1970-01-01", "GMT")),"%d/%m/%Y" ))
# fill NAs in non returning duration
#t_qual$duration <-as.double(round(as.POSIXct(strptime(t_qual$return, "%d/%m/%Y %H:%M"), "GMT") -
#                          as.POSIXct(strptime(t_qual$departure, "%d/%m/%Y %H:%M"), "GMT"), 1))
# write out and update
#write.csv(t_qual, 'C:/seabirds/data/tracking_trip_decisions.csv', quote=F, row.names=F)

# REMOVE 3 unwanted datasets
t_qual<-t_qual[-which(t_qual$ID%in%c('MEND2_RFBO_Christmas', 'MEND3_RFBO_Christmas',
                                     'AUST3_RFBO_Little Cayman', 'GILM10_RFBO_Isabel')),]

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

p2$Longitude<-st_coordinates(p2)[,1]
p2$Latitude<-st_coordinates(p2)[,2]
# **
st_geometry(p2)<-NULL
write.csv(p2, 'C:/seabirds/sourced_data/tracking_data/tracking_master_forage.csv', quote=F, row.names=F)


#calculate scale of interaction using BL scripts

# need col locs
colz<-st_read('C:/seabirds/data/GIS/trackingID_colony_locs.shp')

# make WTSH l/S trip split
t_qual$WTSH_SL<-ifelse(t_qual$duration>(24*3) | t_qual$max_dist>300, 'L', 'S')
lkup<-which(t_qual$sp=='WTSH')
lkup<-lkup[!lkup %in% grep('Aride', t_qual$ID)]
t_qual$ID<-as.character(t_qual$ID)
t_qual[lkup,]$ID<-paste0(t_qual[lkup,]$ID, t_qual[lkup,]$WTSH_SL)

hvals_out<-NULL
for( i in unique(t_qual$ID))#
{
  j<-i
  
  if(i %in% unique(t_qual[lkup,]$ID)){
    i<-substr(i, 1, (nchar(i)-1))}

  cleand<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/clean/',i, '.csv'))
  # filter out bad trips
  cleand<-cleand[cleand$trip_id%in% t_qual[t_qual$ID==j & t_qual$manual_keep=='Y',]$trip_id,]
  cleand$DateTime <- as.POSIXct(strptime(as.character(cleand$DateTime),format="%Y-%m-%d %H:%M:%S",tz="GMT"))
  cleand$ID<-factor(cleand$trip_id) # run per trip rather than track

  s1<-t_qual[t_qual$ID==j,]
  coly<-data.frame(Longitude=as(colz[colz$ID==i,], 'Spatial')@coords[1],
                   Latitude= as(colz[colz$ID==i,], 'Spatial')@coords[2])
  # make into spdf
  sptz <- SpatialPoints(cleand[,4:3],
                        proj4string=CRS("+proj=longlat + datum=wgs84"))
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", coly$Longitude,
                        " +lat_0=", coly$Latitude, sep=""))
  sptz <- spTransform(sptz, CRS=proj.UTM)
  spdf<-SpatialPointsDataFrame(sptz, data=cleand)
  
  if(i== "ZAMO1_BRBO_Pajarera"){next} #not working for this one
  
  # setting res to 2km to match approx grid resolution for kde
  Hvals<-findScale(spdf, scaleARS = T, sumTrips = s1,  res = 2, peakMethod = 'first')
  Hvals$scaleARS_max<-findScale(spdf, scaleARS = T,  res = 2, peakMethod = 'max')$scaleARS
  Hvals$scaleARS_steep<-findScale(spdf, scaleARS = T,  res = 2, peakMethod = 'steep')$scaleARS
  
      # setting res to 2km to match approx grid resolution for kde
      #Hvals<-try(findScale(spdf, scaleARS = T, sumTrips = s1,  res = 2, peakMethod = 'first'))
      #if (class(Hvals)=="try-error"){
      #  Hvals<-data.frame(med_max_dist=NA, mag=NA, scaled_mag=NA, href=NA, ARSscale=NA)}

  hvals_out<-rbind(hvals_out, data.frame(ID=j, Hvals))  
  print(hvals_out)
}

hvals_out$sp<-do.call(c, lapply(strsplit(as.character(hvals_out$ID), '_'), function(x)x[2]))
# lookup interpolation interval
int_sum<-read.csv('C:/seabirds/data/tracking_interval_decisions.csv')
hvals_out<-left_join(hvals_out, int_sum[,c(1, 9)], by='ID')
hvals_out[88:92,]$int_decision<-600
hvals_out[101:102,]$int_decision<-1200
hvals_out[93:100,]$int_decision<-900

hvals_out$sp_group<-hvals_out$sp

hvals_out[hvals_out$sp =='WTSH',]$sp_group<-'WTST' # Aride and all short trips
hvals_out[hvals_out$sp =='WTSH' & substr(hvals_out$ID, nchar(hvals_out$ID), 
          nchar(hvals_out$ID))=='L',]$sp_group<-'WTLG'


hvals_out[hvals_out$sp %in% c('GRFR', 'LEFR', 'MAFR'),]$sp_group<-'FRBD'
hvals_out[hvals_out$sp %in% c('RBTB', 'RTTB'),]$sp_group<-'TRBD'
hvals_out[hvals_out$sp %in% c('BRNO', 'LENO', 'BLNO'),]$sp_group<-'NODD'
hvals_out[hvals_out$sp %in% c('CRTE', 'ROTE', 'CATE'),]$sp_group<-'TERN'

ggplot(data=hvals_out, aes(x=sp_group, y=ARSscale))+
  geom_boxplot()+geom_point(aes(colour=factor(int_decision)), shape=1)

hvals_out%>%group_by(sp_group)%>%summarise(med_hval=median(mag, na.rm=T))
  
write.csv(hvals_out, 'C:/seabirds/data/dataID_hvals.csv', quote=F, row.names=F)


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