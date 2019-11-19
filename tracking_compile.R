# 18/11/19
# Code to compile contributed tracking datasets into standardised master csv
# start with .csv files and then include more complex

# WTSH data from Congdon lab - set up master
#congdon_wtsh<-read.csv('C:/seabirds/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm.csv')
#master<-data.frame(
#dataID=rep('CONG1', nrow(congdon_wtsh)),
#sp=rep('WTSH', nrow(congdon_wtsh)),
#colony=congdon_wtsh$Colony,
#trackID=congdon_wtsh$TrackID,
#date=congdon_wtsh$Date,
#time=congdon_wtsh$Time,
#latitude=congdon_wtsh$Latitude,
#longitude=congdon_wtsh$Longitude)
#write.csv(master, 'C:/seabirds/sourced_data/tracking_data/tracking_master.csv', quote=F, row.names=F)

library(sf)
master<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master.csv')

# Add Congdon MABO and BRBO data
congdon_brbo<-read.csv('C:/seabirds/phd/analyses/BRBO_raine/BRBO_raine_hmm.csv')
master<-rbind(master, data.frame(dataID='CONG2',sp='BRBO',colony='Raine',
  trackID=congdon_brbo$TrackID, date=congdon_brbo$Date,time=congdon_brbo$Time,
  latitude=congdon_brbo$Latitude,longitude=congdon_brbo$Longitude))

congdon_mabo<-read.csv('C:/coral_fish/teaching/supervision/Charlotte_Dean/data/Swains_2014_tracking_trips_clean.csv')
master<-rbind(master, data.frame(dataID='CONG3',sp=congdon_mabo$SpeciesID, colony='Swains',
              trackID=congdon_mabo$TrackID,date=congdon_mabo$Date,time=congdon_mabo$Time,
              latitude=congdon_mabo$Latitude,longitude=congdon_mabo$Longitude))

congdon_mabo<-read.csv('C:/coral_fish/teaching/supervision/Charlotte_Dean/data/Swains_2015_tracking_trips_clean.csv')
master<-rbind(master, data.frame(dataID='CONG4',sp=congdon_mabo$SpeciesID, colony='Swains',
              trackID=congdon_mabo$TrackID, date=congdon_mabo$Date,time=congdon_mabo$Time,
              latitude=congdon_mabo$Latitude,longitude=congdon_mabo$Longitude))

# now pull some stranger formatted datasets

# Surman Lesser Noddy
p1<-st_read('C:/seabirds/sourced_data/tracking_data/raw/Surman_WA_LSNO/all_tracks.shp')
p1<-as.data.frame(p1)
p2<-strsplit(paste(p1$geometry), ',')
p2_out<-NULL
for(i in 1:length(p2))
  {p2_out<-rbind(p2_out, data.frame(latitude=as.numeric(p2[[i]][2]),
                                    longitude=substr(p2[[i]][1], 3, nchar(p2[[i]][1]))))}

master<-rbind(master, data.frame(dataID='SURM1',sp='LENO', colony='Pelsaert',
                                 trackID=p1$layer, date=substr(p1$timestamp, 1,10),
                                 time=substr(p1$timestamp, 12,19),p2_out))

