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
#date=gsub('-', '/', substr(congdon_wtsh$DateTime, 1, 10)),
#time=substr(congdon_wtsh$DateTime, 12, 19),
#latitude=congdon_wtsh$Latitude,
#longitude=congdon_wtsh$Longitude)
#write.csv(master, 'C:/seabirds/sourced_data/tracking_data/tracking_master.csv', quote=F, row.names=F)

library(sf)
library(readxl)
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

# now pull the other datasets

# Surman Lesser Noddy
p1<-st_read('C:/seabirds/sourced_data/tracking_data/raw/Surman_WA_LSNO/all_tracks.shp')
p1<-as.data.frame(p1)
p2<-strsplit(paste(p1$geometry), ',')
p2_out<-NULL
for(i in 1:length(p2))
  {p2_out<-rbind(p2_out, data.frame(latitude=as.numeric(p2[[i]][2]),
                                    longitude=as.numeric(substr(p2[[i]][1], 3, nchar(p2[[i]][1])))))}

master<-rbind(master, data.frame(dataID='SURM1',sp='LENO', colony='Pelsaert',
                                 trackID=p1$layer, date=substr(p1$timestamp, 1,10),
                                 time=substr(p1$timestamp, 12,19),p2_out))

# Neumann Sooty Terns
p1<-NULL
for( i in 2:7)# exclude 1 as 6 lines
{p1<-rbind(p1, read_xlsx('C:/seabirds/sourced_data/tracking_data/raw/Neumann_Seychelles_SOTE.xlsx', sheet=i))}

master<-rbind(master, data.frame(dataID='NEUM1',sp='SOTE', colony='Seychelles',
                                 trackID=p1$Tag,date=gsub('-', '/', p1$Date),time=p1$Time_T,
                                 latitude=p1$Y,longitude=p1$X))

# Shephard Brown Noddy and Sooty Tern

p1<-read_xlsx('C:/seabirds/sourced_data/tracking_data/raw/Shephard_WA_BRNO.xlsx', sheet=3)
p1<-p1[1:12812,]
p1$sp='BRNO'
p1[which(p1$Species==2),]$sp<-'SOTE'
p1$colony='Lancelin'
p1[which(p1$Species==2),]$colony<-'Rat'
master<-rbind(master, data.frame(dataID='SHEP1',sp=p1$sp, colony=p1$colony,
                                 trackID=factor(p1$ID),date=gsub('-', '/', substr(p1$DateTimeWST, 1,10)),time=substr(p1$DateTimeWST, 12, 19),
                                 latitude=p1$LatDD,longitude=p1$LongDD))

# Yoda BRBO
birds<-list.files('C:/seabirds/sourced_data/tracking_data/raw/BrownBoobyGPS_Yoda')
birds2<-NULL
for(i in 2:length(birds)) # skip logger 1 as different format
{birds2<-rbind(birds2, read.csv(paste0('C:/seabirds/sourced_data/tracking_data/raw/BrownBoobyGPS_Yoda/', birds[i])))}

#sprintf handy formatting
birds2<-data.frame(dataID='YODA1',sp='BRBO', colony='Nakanokamishima',
                   trackID=birds2$Logger.ID,date=paste(birds2$Year, sprintf("%02d", birds2$Month), sprintf("%02d", birds2$Day), sep='/'),
                   time=paste(sprintf("%02d",birds2$Hour), sprintf("%02d", birds2$Minute), sprintf("%02d", birds2$Second), sep=':'),
                   latitude=birds2$Latitude,longitude=birds2$Longitude)

bird1<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/raw/BrownBoobyGPS_Yoda/', birds[1]))
birds2<-rbind(birds2, data.frame(dataID='YODA1',sp='BRBO', colony='Nakanokamishima',
           trackID='GiPSy',date=substr(bird1$Date.Time, 1, 10),
           time=substr(bird1$Date.Time, 12, 19),
           latitude=bird1$Latitude,longitude=bird1$Longitude))
master<-rbind(master, birds2)

# Bunce BRBO
b1<-read_xls('C:/seabirds/sourced_data/tracking_data/raw/Bunce_Swains_BRBO/G53_060125085537.xls', sheet='all data')
b1$trackID='G53_060125085537'
b2<-read_xls('C:/seabirds/sourced_data/tracking_data/raw/Bunce_Swains_BRBO/G54_060125090457.xls', sheet='all data')
b2$trackID='G54_060125090457'
b3<-read_xls('C:/seabirds/sourced_data/tracking_data/raw/Bunce_Swains_BRBO/G55_060125090823.xls', sheet='all data')
b3$trackID='G55_060125090823'
b4<-read_xls('C:/seabirds/sourced_data/tracking_data/raw/Bunce_Swains_BRBO/G102_060125085121B.xls', sheet='all data')
b4$trackID='G102_060125085121B'
b5<-read_xls('C:/seabirds/sourced_data/tracking_data/raw/Bunce_Swains_BRBO/G104_060125084110B.xls', sheet='all data')
b5$trackID='G104_060125084110B'
b1<-rbind(b1, b2, b3)
b2<-rbind(b4, b5)
b3<-rbind(b1[,c(13, 1, 2, 6, 10)], b2[,c(10, 1, 2, 4, 6)])

master<-rbind(master, data.frame(dataID='BUNC1',sp='BRBO', colony='Swains',
                                 trackID=b3$trackID,date=paste(substr(b3$DATE, 7,10), substr(b3$DATE, 4,5), substr(b3$DATE, 1,2), sep='/'),
                                 time=substr(b3$TIME, 12, 19),
                                 latitude=b3$Y,longitude=b3$X))

# Machovsy-Capuska MABO
birds<-list.files('C:/seabirds/sourced_data/tracking_data/raw/Capuska_MABO', recursive=T) # recursive NICE!
birds2<-NULL
for(i in 1:length(birds))
{birds2<-rbind(birds2, read.csv(paste0('C:/seabirds/sourced_data/tracking_data/raw/Capuska_MABO/', birds[i]), h=F))}

# Note there are some 0,0 lat longs in data
master<-rbind(master, data.frame(dataID='CAPU1',sp='MABO', colony='LHI',
                                 trackID=factor(birds2$V2),
                                 date=paste('2013', unlist(lapply(strsplit(as.character(birds2$V3), '\\.'),
                                                    function(x){sprintf("%02d",as.numeric(x[2]))})),
                                                    unlist(lapply(strsplit(as.character(birds2$V3), '\\.'),
                                                    function(x){sprintf("%02d",as.numeric(x[1]))})), sep='/'),
                                 time=birds2$V16, latitude=birds2$V7,longitude=birds2$V6))

# Ravache WTSH
p1<-read.table('C:/seabirds/sourced_data/tracking_data/raw/Ravache_Newcal_Ardenna_Pacifica_Trips_Caledonia_171819.txt', sep='\t', h=T)

master<-rbind(master, data.frame(dataID='RAVA1',sp='WTSH', colony=p1$Site,
                                 trackID=p1$ID,date=paste(substr(p1$TIME_LOC, 7,10), substr(p1$TIME_LOC, 4,5), substr(p1$TIME_LOC, 1,2), sep='/'),
                                 time=paste0(substr(p1$TIME_LOC, 12, 16),':00'),
                                 latitude=p1$Latitude,longitude=p1$Longitude))

