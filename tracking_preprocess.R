# tracking clean
library(dplyr)
library(ggplot2)
library(adehabitatLT)
# soure BL scripts
source('C:/seabirds/sourced_data/code/new_mIBA/formatFields.R')
source('C:/seabirds/sourced_data/code/new_mIBA/tripSplit.R')


# make summary table of tracking interval per dataset

master<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master.csv')
# Capuska data 0 lat long  quick fix
master<-master[-which(master$latitude==0 & master$longitude==0),]

master$ID<-paste(master$dataID, master$sp, master$colony, sep='_')

#master<-master%>%group_by(ID, trackID)%>%mutate(dt=c(0, diff(tracktime)))

#master_sum<-master%>%group_by(ID)%>%summarise(npoints=n(), nBird=length(unique(trackID)),
#                                              med_int=median(dt),
#                                              int_25=quantile(dt, 0.25, na.rm=T),
#                                              int_75=quantile(dt,0.75, na.rm=T))
#write.csv(master_sum, 'C:/seabirds/data/tracking_interval_decisions.csv', quote=F, row.names=F)
#p1<-ggplot(data=master, aes(x=dt))+geom_histogram()+facet_wrap(~ID, scales='free')
#png('C:/seabirds/results/tracking_intervals.png',width = 12, height =12 , units ="in", res =600)
#p1
#dev.off()

# Data processing
# loop to:
# remove speeds > 75km/h
# interpolate to set interval
# perform tripsplit
 
# read in interval assigned data summary
int_sum<-read.csv('C:/seabirds/data/tracking_interval_decisions.csv')

for(i in int_sum$ID)
{
  #select dataset
  d1<-master[master$ID==i,]
  # fix for dataID=='OPPE1'
  #if(i=='OPPE1_MABO_St Helena'){
  # d1[d1$dt<2300000,]$dt<-NA 
  # d1<-tidyr::fill(d1, dt, .direction='down')
  # d1$trackID<-paste(d1$trackID, as.character(d1$dt))
  #}
  
  # recalc tracktime just to make sure
  d1$datetime <- as.POSIXct(strptime(paste(d1$date, d1$time, sep=""), "%Y/%m/%d %H:%M:%S"), "GMT")
  d1$tracktime <- as.double(d1$datetime)
  
  # get approx colony loc
  trk_start<-d1%>%group_by(trackID)%>%summarise_all(first)
  # project data to lambert equal area centred on data
  ### CREATE PROJECTED DATAFRAME ###
  sptz <- SpatialPoints(data.frame(d1$longitude, d1$latitude),
                             proj4string=CRS("+proj=longlat + datum=wgs84"))
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", median(trk_start$longitude),
                        " +lat_0=", median(trk_start$latitude), sep=""))
  sptz <- spTransform(sptz, CRS=proj.UTM)
  
  d1$trackID<-as.character(d1$trackID) 
  if("" %in% unique(d1$trackID)){d1[d1$trackID=="",]$trackID<-'noid'}
  
  trajectories <- as.ltraj(xy=data.frame(sptz@coords[,1],
                                         sptz@coords[,2]),
                           date=as.POSIXct(d1$tracktime, origin="1970/01/01", tz="GMT"),
                           id=d1$trackID, typeII = TRUE)   
  
  rm(sptz)#save mem
  #to data.frame
  t1<-ld(trajectories)
  if(length(which((t1$dist/t1$dt)>25))>0){
  t1<-t1[-which((t1$dist/t1$dt)>25),]} #~ 90 km/h Mendez
  # back to ltraj
  trajectories<-dl(t1)
  #interpolate to lookup val
  trajectories <- redisltraj(trajectories, int_sum[int_sum$ID==i,]$int_decision, type="time")
  
  # to data.frame
  trajectories<-ld(trajectories)
  #points back to latlong
  sptz <- SpatialPoints(data.frame(trajectories$x, trajectories$y),
                        proj4string=proj.UTM)
  sptz <- spTransform(sptz, CRS=CRS("+proj=longlat + datum=wgs84"))
  trajectories<-data.frame(trajectories, sptz@coords)
  
  trajectories<-formatFields(trajectories, field_ID   = "id", field_DateTime='date',
                     field_Lon  = "trajectories.x", field_Lat  = "trajectories.y")
 #colony
  coly<-data.frame(Longitude=median(trk_start$longitude),Latitude= median(trk_start$latitude))
  
  #tripsplit
  dat_int_trips<-tripSplit(tracks = trajectories, Colony= coly, 
                        InnerBuff  = 5,ReturnBuff = 20, Duration   = 1,
                        plotit     = F,  rmColLocs  = T, cleanDF = T)
  
  #write out with interpolation
  write.csv(dat_int_trips@data, paste0('C:/seabirds/sourced_data/tracking_data/clean/',
                                      i, '.csv'), quote=F, row.names=F)

  rm(dat_int_trips)#save mem
  print(i)
  print(paste(which(i==int_sum$ID),'of', length(int_sum$ID)))
}

# LOOP to tripsplit first! - only ran for for OPPE1_MABO_St Helena

# read in interval assigned data summary
int_sum<-read.csv('C:/seabirds/data/tracking_interval_decisions.csv')

for(i in int_sum$ID)
{
  #select dataset
  d1<-master[master$ID==i,]
  
  # recalc tracktime just to make sure
  d1$datetime <- as.POSIXct(strptime(paste(d1$date, d1$time, sep=""), "%Y/%m/%d %H:%M:%S"), "GMT")
  d1$tracktime <- as.double(d1$datetime)
  
  # get approx colony loc
  trk_start<-d1%>%group_by(trackID)%>%summarise_all(first)
  
  #colony
  coly<-data.frame(Longitude=median(trk_start$longitude),Latitude= median(trk_start$latitude))
  
  d1$ID<-NULL # rm for formatFields()
  d1<-formatFields(d1, field_ID   = "trackID", field_DateTime='datetime',
                             field_Lon  = "longitude", field_Lat  = "latitude")

  #tripsplit
  dat_int_trips<-tripSplit(tracks = d1, Colony= coly, 
                           InnerBuff  = 3,ReturnBuff = 10, Duration   = 1,
                           plotit     = F,  rmColLocs  = T, cleanDF = T)
  
  
  dat_int_trips$trip_id<-as.character(dat_int_trips$trip_id) 
  if("" %in% unique(dat_int_trips$trip_id)){dat_int_trips[dat_int_trips$trip_id=="",]$trip_id<-'noid'}
  
  trajectories <- as.ltraj(xy=data.frame(dat_int_trips@coords[,1],
                                         dat_int_trips@coords[,2]),
                           date=dat_int_trips$DateTime,
                           id=dat_int_trips$trip_id, typeII = TRUE)   
  
 
  #to data.frame
  t1<-ld(trajectories)
  if(length(which((t1$dist/t1$dt)>25))>0){
    print(paste('removed', length(which((t1$dist/t1$dt)>25)),
                'over 90kmph of total', nrow(t1)))
    t1<-t1[-which((t1$dist/t1$dt)>25),]} #~ 90 km/h Mendez
  # back to ltraj
  trajectories<-dl(t1)
  #interpolate to lookup val
  trajectories <- redisltraj(trajectories, int_sum[int_sum$ID==i,]$int_decision, type="time")
  
  # to data.frame
  trajectories<-ld(trajectories)
  #points back to latlong
  sptz <- SpatialPoints(data.frame(trajectories$x, trajectories$y),
                        proj4string=dat_int_trips@proj4string)
  sptz <- spTransform(sptz, CRS=CRS("+proj=longlat + datum=wgs84"))
  trajectories<-data.frame(trajectories, sptz@coords)
  
 
  #write out with interpolation
  write.csv(trajectories, paste0('C:/seabirds/sourced_data/tracking_data/clean_alt/',
                                       i, '.csv'), quote=F, row.names=F)
  
  rm(dat_int_trips)#save mem
  print(i)
  print(paste(which(i==int_sum$ID),'of', length(int_sum$ID)))
}
