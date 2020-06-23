# tracking clean
library(dplyr)
library(ggplot2)
library(adehabitatLT)
library(sf)
library(maptools)
# soure BL scripts
source('C:/seabirds/sourced_data/code/new_mIBA/formatFields.R')
source('C:/seabirds/sourced_data/code/new_mIBA/tripSplit.R')
source('C:/seabirds/sourced_data/code/new_mIBA/tripSummary.R')


# make summary table of tracking interval per dataset

master<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master.csv')

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

# get colony locs for each ID, file will be edited in QGIS
#sf_col<-master%>%group_by(ID, trackID)%>%summarise_all(first)%>%ungroup()%>%
#group_by(ID)%>%summarise(latitude=median(latitude), longitude=median(longitude))%>%
#  st_as_sf(coords=c('longitude', 'latitude'), crs=4326)
#st_write(sf_col, 'C:/seabirds/data/GIS/trackingID_colony_locs.shp')

# to assist the first and last 5 points of each track
#sf_trax<-master%>%group_by(ID, trackID)%>%slice(1:10)%>%
#  dplyr::select(latitude, longitude)%>%
#  st_as_sf(coords=c('longitude', 'latitude'), crs=4326)
#st_write(sf_trax, 'C:/seabirds/data/GIS/tracking_ID_firstpoints.shp')

# Data processing
# loop to:
# remove speeds > 90km/h
# interpolate to set interval
# perform tripsplit
 
# read in interval assigned data summary
int_sum<-read.csv('C:/seabirds/data/tracking_interval_decisions.csv')
# read in col locs
colz<-st_read('C:/seabirds/data/GIS/trackingID_colony_locs.shp')

for(i in int_sum$ID)
{
  #select dataset
  d1<-master[master$ID==i,]
   # edited dataset with near colony loitering removed

  if(i=='OPPE1_MABO_St Helena'){next} # needs to be tripsplit first

  # recalc tracktime just to make sure
  d1$datetime <- as.POSIXct(strptime(paste(d1$date, d1$time, sep=""), "%Y/%m/%d %H:%M:%S"), "GMT")
  d1$tracktime <- as.double(d1$datetime)
  
  #colony from lookup
  coly<-data.frame(Longitude=as(colz[colz$ID==i,], 'Spatial')@coords[1],
                   Latitude= as(colz[colz$ID==i,], 'Spatial')@coords[2])
  
  # project data to lambert equal area centred on data
  ### CREATE PROJECTED DATAFRAME ###
  sptz <- SpatialPoints(data.frame(d1$longitude, d1$latitude),
                             proj4string=CRS("+proj=longlat + datum=wgs84"))
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", coly$Longitude,
                        " +lat_0=", coly$Latitude, sep=""))
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
  print(paste('nrow pre speed filt', nrow(t1)))
  if(length(which((t1$dist/t1$dt)>25))>0){
  t1<-t1[-which((t1$dist/t1$dt)>25),]} #~ 90 km/h Mendez
  if(i=='AUST4_MAFR_Little Cayman'){t1<-ld(trajectories)}# fix to not remove outlers from this dset
  print(paste('nrow post speed filt', nrow(t1)))
  # back to ltraj
  trajectories<-dl(t1)
  #interpolate to lookup val
  # Don't interpolate this dataset, gives better embc results
  if(i=='AUST4_MAFR_Little Cayman'){}else{
  trajectories <- redisltraj(trajectories, int_sum[int_sum$ID==i,]$int_decision, type="time")}
  
  # to data.frame
  trajectories<-ld(trajectories)
  #points back to latlong
  sptz <- SpatialPoints(data.frame(trajectories$x, trajectories$y),
                        proj4string=proj.UTM)
  sptz <- spTransform(sptz, CRS=CRS("+proj=longlat + datum=wgs84"))
  trajectories<-data.frame(trajectories, sptz@coords)
  
  trajectories<-formatFields(trajectories, field_ID   = "id", field_DateTime='date',
                     field_Lon  = "trajectories.x", field_Lat  = "trajectories.y")

  #tripsplit
  dat_int_trips<-tripSplit(tracks = trajectories, Colony= coly, 
                        InnerBuff  = 3,ReturnBuff = 10, Duration   = 1,
                        plotit     = F,  rmColLocs  = T, cleanDF = T)
  
  # fix for multi colony
  if(i=='AUST1_BRBO_Cayman Brac'){
    dif_tripz<-c('B01', 'B02', 'B06', 'B12', 'B13', 'C05', 'C25','C28', 'P00', 'P03', 'P04', 'P07', 'P08')
    dok<-dat_int_trips[-which(dat_int_trips$ID %in% dif_tripz),]
    tempy1<-tripSplit(tracks = trajectories[trajectories$ID %in% dif_tripz,],
                      Colony= data.frame(Longitude=-79.7293, Latitude= 19.7549), 
                      InnerBuff  = 3,ReturnBuff = 10, Duration   = 1,
                      plotit     = F,  rmColLocs  = T, cleanDF = T)
     tempy1@proj4string<-dat_int_trips@proj4string
       dat_int_trips<-spRbind(dok,tempy1)}
  
  #write out with interpolation
  write.csv(dat_int_trips@data, paste0('C:/seabirds/sourced_data/tracking_data/clean/',
                                      i, '.csv'), quote=F, row.names=F)

  rm(dat_int_trips)#save mem
  print(i)
  print(paste(which(i==int_sum$ID),'of', length(int_sum$ID)))
}

## Loop to assist with the manual removal of crap trips
## Loop to read raw and clean datasets produce per tripsplit figure

# read in interval assigned data summary
int_sum<-read.csv('C:/seabirds/data/tracking_interval_decisions.csv')
# read in col locs
colz<-st_read('C:/seabirds/data/GIS/trackingID_colony_locs.shp')

trip_table<-NULL
for(i in int_sum$ID)
{
  #select dataset
  rawd<-master[master$ID==i,]
  
  #write out raw
  write.csv(rawd, paste0('C:/seabirds/sourced_data/tracking_data/compiled/',
                                       i, '.csv'), quote=F, row.names=F)
  
  if(i=='OPPE1_MABO_St Helena'){next} # different columns
  
  cleand<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/clean/',i, '.csv'))
  
  cleand <- SpatialPointsDataFrame(SpatialPoints(data.frame(cleand$Longitude, 
            cleand$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84")),
            data = cleand, match.ID=F)
  mid_point<-data.frame(geosphere::centroid(cbind(cleand$Longitude, cleand$Latitude)))
  
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
  cleand <- spTransform(cleand, CRS=proj.UTM)
  cleand$DateTime <- as.POSIXct(strptime(cleand$DateTime, "%Y-%m-%d %H:%M:%S"), "GMT")
  
  coly<-data.frame(Longitude=as(colz[colz$ID==i,], 'Spatial')@coords[1],
                   Latitude= as(colz[colz$ID==i,], 'Spatial')@coords[2])
  
  #tripSummary
  sumTrips <- tripSummary(Trips = cleand, Colony = coly)
  sumTrips<-rename(sumTrips, trackID=ID)
  
  # augment, rbind and writeout summary
  sum_trips<-data.frame(ID=i, sumTrips)

  trip_table<-rbind(trip_table, sum_trips)
  write.csv(trip_table, 'C:/seabirds/data/tracking_trip_decisions.csv', quote=F, row.names=F)
  
   # output plots
  
  trackz<-unique(cleand$ID)
  
  if(length(trackz)<25){
  
  TRACKPLOT <- cleand@data %>%
    arrange(.data$ID, .data$TrackTime) %>%
    ggplot(aes(., x=Longitude, y=Latitude, colour=substr(trip_id,
                  nchar(as.character(trip_id)), nchar(as.character(trip_id))))) +
    geom_path(colour='grey') +
    geom_point(size=0.2)+
    geom_point(data=coly, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
    facet_wrap(ggplot2::vars(ID), scales='free') +
    theme(panel.background=element_rect(fill="white", colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour="black", fill="white"),
          panel.border = element_blank(), legend.position = "none")
  
    png(paste0('C:/seabirds/sourced_data/tracking_data/fig_clean_tripsplit/',i,'.png'),width = 12, height =12 , units ="in", res =600)
    print(TRACKPLOT)
    dev.off()
    }else{lapply(split(trackz, ceiling(seq_along(trackz)/25)),
                      function(x){
                                TRACKPLOT <- cleand@data %>% filter(ID %in% x)%>%
                                  arrange(.data$ID, .data$TrackTime) %>%
                                  ggplot(aes(., x=Longitude, y=Latitude, colour=substr(trip_id,
                                                                                       nchar(as.character(trip_id)), nchar(as.character(trip_id))))) +
                                  geom_path(colour='grey') +
                                  geom_point(size=0.2)+
                                  geom_point(data=coly, aes(x=Longitude, y=Latitude), col='red', shape=16, size=2) +
                                  facet_wrap(ggplot2::vars(ID), scales='free') +
                                  theme(panel.background=element_rect(fill="white", colour="black"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        strip.background = element_rect(colour="black", fill="white"),
                                        panel.border = element_blank(), legend.position = "none")
                                
                                png(paste0('C:/seabirds/sourced_data/tracking_data/fig_clean_tripsplit/',i,
                                           '_', x[[1]][1], '.png'),width = 12, height =12 , units ="in", res =600)
                                print(TRACKPLOT)
                                dev.off()})
                                } # end of else loop
 
print(i)
print(paste(which(i==int_sum$ID),'of', length(int_sum$ID)))      
} # end of big loop
  
# add breedstage
t2<-left_join(trip_table, master%>%group_by(trackID)%>%summarise_all(first)%>%dplyr::select(trackID,breedstage), by='trackID')
write.csv(t2, 'C:/seabirds/data/tracking_trip_decisions.csv', quote=F, row.names=F)


# Section to resample some datasets to match others once we collapse by colony name

for( i in c ('MCDU1_MABO_Swains', 'MCDU1_BRBO_Swains', 'CONG3_BRBO_Swains', 'SOAN6_BRBO_Dog'))
{
  # note original data stored in clean_alt folder
  cleand<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/clean_alt/',i, '.csv'))
  
  #colony from lookup
  coly<-data.frame(Longitude=as(colz[colz$ID==i,], 'Spatial')@coords[1],
                   Latitude= as(colz[colz$ID==i,], 'Spatial')@coords[2])
  
  sptz <- SpatialPoints(data.frame(cleand$Longitude, cleand$Latitude),
                        proj4string=CRS("+proj=longlat + datum=wgs84"))
  proj.UTM <- CRS(paste("+proj=laea +lon_0=", coly$Longitude,
                        " +lat_0=", coly$Latitude, sep=""))
  sptz <- spTransform(sptz, CRS=proj.UTM)
  
  trajectories <- as.ltraj(xy=data.frame(sptz@coords[,1],
                                         sptz@coords[,2]),
                           date=as.POSIXct(cleand$TrackTime, origin="1970/01/01", tz="GMT"),
                           id=cleand$trip_id, typeII = TRUE)   
  
  if(i=='MCDU1_MABO_Swains'){inty=60}
  if(i=='MCDU1_BRBO_Swains'){inty=120}
  if(i=='CONG3_BRBO_Swains'){inty=120}
  if(i=='SOAN6_BRBO_Dog'){inty=120}
  trajectories <- redisltraj(trajectories, inty, type="time")
  
  # to data.frame
  trajectories<-ld(trajectories)
  #points back to latlong
  sptz <- SpatialPoints(data.frame(trajectories$x, trajectories$y),
                        proj4string=proj.UTM)
  sptz <- spTransform(sptz, CRS=CRS("+proj=longlat + datum=wgs84"))
  trajectories<-data.frame(trajectories, sptz@coords)
  
  trajectories<-formatFields(trajectories, field_ID   = "id", field_DateTime='date',
                             field_Lon  = "trajectories.x", field_Lat  = "trajectories.y")
  
  #tripsplit
  dat_int_trips<-tripSplit(tracks = trajectories, Colony= coly, 
                           InnerBuff  = 3,ReturnBuff = 10, Duration   = 1,
                           plotit     = F,  rmColLocs  = T, cleanDF = T)
  
  dat_int_trips$ID<-substr(dat_int_trips$ID, 1, (nchar(as.character(dat_int_trips$ID))-1))
  dat_int_trips$trip_id<-substr(dat_int_trips$trip_id, 1, (nchar(as.character(dat_int_trips$trip_id))-1))
  
  #write out with interpolation
  write.csv(dat_int_trips@data, paste0('C:/seabirds/sourced_data/tracking_data/clean/',
                                       i, '.csv'), quote=F, row.names=F)
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
