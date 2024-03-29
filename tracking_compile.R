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
#longitude=congdon_wtsh$Longitude,
#breedstage='chick')
#write.csv(master, 'C:/seabirds/sourced_data/tracking_data/tracking_master.csv', quote=F, row.names=F)

# data req
#wtsh<-filter(congdon_wtsh, Colony=='LHI')
#d1<-data.frame(Common_name='Wedge-tailed Shearwater', Site_name='Lord Howe Island', Colony_name='Signal Point', Lat_colony=-31.52459, Long_colony=159.05991,
#               Device='GPS', Bird_id=wtsh$TrackID, Track_id=wtsh$trip_id, Sex='unknown', Breed_stage='chick', Breed_status='Breeding',
#               Date_AEST=wtsh$DateAEST, Time_AEST=wtsh$TimeAEST, Longitude=wtsh$Longitude, Latitude=wtsh$Latitude, Argos_quality=NA)
#write.csv(d1, 'C:/seabirds/admin/Clay_wind_collab/wtsh_lhi.csv', quote=F, row.names=F )

library(sf)
library(readxl)
library(dplyr)
#master<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master.csv')

# Add Congdon MABO and BRBO data
congdon_brbo<-read.csv('C:/seabirds/phd/analyses/BRBO_raine/BRBO_raine_hmm.csv')
master<-rbind(master, data.frame(dataID='CONG2',sp='BRBO',colony='Raine',
  trackID=congdon_brbo$TrackID, date=congdon_brbo$Date,time=congdon_brbo$Time,
  latitude=congdon_brbo$Latitude,longitude=congdon_brbo$Longitude,breedstage='chick'))

congdon_mabo<-read.csv('C:/coral_fish/teaching/supervision/Charlotte_Dean/data/Swains_2014_tracking_trips_clean.csv')
master<-rbind(master, data.frame(dataID='CONG3',sp=congdon_mabo$SpeciesID, colony='Swains',
              trackID=congdon_mabo$TrackID,date=congdon_mabo$Date,time=trimws(congdon_mabo$Time),
              latitude=congdon_mabo$Latitude,longitude=congdon_mabo$Longitude,breedstage='chick'))

congdon_mabo<-read.csv('C:/coral_fish/teaching/supervision/Charlotte_Dean/data/Swains_2015_tracking_trips_clean.csv')
master<-rbind(master, data.frame(dataID='CONG4',sp=congdon_mabo$SpeciesID, colony='Swains',
              trackID=congdon_mabo$TrackID, date=congdon_mabo$Date,time=congdon_mabo$Time,
              latitude=congdon_mabo$Latitude,longitude=congdon_mabo$Longitude,breedstage='chick'))

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
                                 time=substr(p1$timestamp, 12,19),p2_out, breedstage='chick'))

# Neumann Sooty Terns
p1<-NULL
for( i in 2:7)# exclude 1 as 6 lines
{p1<-rbind(p1, read_xlsx('C:/seabirds/sourced_data/tracking_data/raw/Neumann_Seychelles_SOTE.xlsx', sheet=i))}

master<-rbind(master, data.frame(dataID='NEUM1',sp='SOTE', colony='Seychelles',
                                 trackID=p1$Tag,date=gsub('-', '/', p1$Date),time=p1$Time_T,
                                 latitude=p1$Y,longitude=p1$X, breedstage='incubation'))

# Shephard Brown Noddy and Sooty Tern

p1<-read_xlsx('C:/seabirds/sourced_data/tracking_data/raw/Shephard_WA_BRNO.xlsx', sheet=3)
p1<-p1[1:12812,]
p1$sp='BRNO'
p1[which(p1$Species==2),]$sp<-'SOTE'
p1$colony='Lancelin'
p1[which(p1$Species==2),]$colony<-'Rat'
p1$breedstage='incubation'
p1[which(p1$Species==2),]$breedstage<-'chick'
master<-rbind(master, data.frame(dataID='SHEP1',sp=p1$sp, colony=p1$colony,
                                 trackID=factor(p1$ID),date=gsub('-', '/', substr(p1$DateTimeWST, 1,10)),time=substr(p1$DateTimeWST, 12, 19),
                                 latitude=p1$LatDD,longitude=p1$LongDD, breedstage=p1$breedstage))

# Yoda BRBO
birds<-list.files('C:/seabirds/sourced_data/tracking_data/raw/BrownBoobyGPS_Yoda')
birds2<-NULL
for(i in 2:length(birds)) # skip logger 1 as different format
{birds2<-rbind(birds2, read.csv(paste0('C:/seabirds/sourced_data/tracking_data/raw/BrownBoobyGPS_Yoda/', birds[i])))}

#sprintf handy formatting
birds2<-data.frame(dataID='YODA1',sp='BRBO', colony='Nakanokamishima',
                   trackID=birds2$Logger.ID,date=paste(birds2$Year, sprintf("%02d", birds2$Month), sprintf("%02d", birds2$Day), sep='/'),
                   time=paste(sprintf("%02d",birds2$Hour), sprintf("%02d", birds2$Minute), sprintf("%02d", birds2$Second), sep=':'),
                   latitude=birds2$Latitude,longitude=birds2$Longitude,  breedstage='chick')

bird1<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/raw/BrownBoobyGPS_Yoda/', birds[1]))
birds2<-rbind(birds2, data.frame(dataID='YODA1',sp='BRBO', colony='Nakanokamishima',
           trackID='GiPSy',date=substr(bird1$Date.Time, 1, 10),
           time=substr(bird1$Date.Time, 12, 19),
           latitude=bird1$Latitude,longitude=bird1$Longitude,  breedstage='chick'))
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
                                 latitude=b3$Y,longitude=b3$X,  breedstage='chick'))

# Machovsy-Capuska MABO
birds<-list.files('C:/seabirds/sourced_data/tracking_data/raw/Capuska_MABO', recursive=T) # recursive NICE!
birds2<-NULL
for(i in 1:length(birds))
{birds2<-rbind(birds2, read.csv(paste0('C:/seabirds/sourced_data/tracking_data/raw/Capuska_MABO/', birds[i]), h=F))}

# Note there are some 0,0 lat longs in data - remove
birds2<-birds2[birds2$V10=='A',]
#remove junk trips
birds2<-birds2[birds2$V17 %in% names(which(table(birds2$V17)>100)),]

master<-rbind(master, data.frame(dataID='CAPU1',sp='MABO', colony='LHI',
                                 trackID=factor(paste(birds2$V17, birds2$V18)),
                                 date=paste('2013', unlist(lapply(strsplit(as.character(birds2$V14), '\\.'),
                                                    function(x){sprintf("%02d",as.numeric(x[2]))})),
                                                    unlist(lapply(strsplit(as.character(birds2$V14), '\\.'),
                                                    function(x){sprintf("%02d",as.numeric(x[1]))})), sep='/'),
                                 time=birds2$V16, latitude=birds2$V7,longitude=birds2$V6,  breedstage='chick'))

# Ravache WTSH
p1<-read.table('C:/seabirds/sourced_data/tracking_data/raw/Ravache_Newcal_Ardenna_Pacifica_Trips_Caledonia_171819.txt', sep='\t', h=T)

master<-rbind(master, data.frame(dataID='RAVA1',sp='WTSH', colony=p1$Site,
                                 trackID=p1$ID,date=paste(substr(p1$TIME_LOC, 7,10), substr(p1$TIME_LOC, 4,5), substr(p1$TIME_LOC, 1,2), sep='/'),
                                 time=paste0(substr(p1$TIME_LOC, 12, 16),':00'),
                                 latitude=p1$Latitude,longitude=p1$Longitude,  breedstage='chick'))
# Zamora BRBO & RBTB

p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Zamora_Mexico_BRBO_RBTB.csv', h=T)
zam_meta<-read_xlsx('C:/seabirds/sourced_data/tracking_data/raw/META_Zamora_Mexico_BRBO_RBTB.xlsx', sheet=1)

p1$sp<-'BRBO'
p1[which(p1$Species=='Phaethon athereus'),]$sp<-'RBTB'
p1$breedstage='incubation'
p1[which(p1$ID %in% zam_meta[which(zam_meta$`Brood size`==1),]$ID),]$breedstage<-'chick'

p1$Date<-as.character(p1$Date)
p1[which(p1$Species=='Phaethon athereus'),]$Date<-paste(substr(p1[which(p1$Species=='Phaethon athereus'),]$Date, 7,10), 
                                                        substr(p1[which(p1$Species=='Phaethon athereus'),]$Date, 4,5), 
                                                        substr(p1[which(p1$Species=='Phaethon athereus'),]$Date, 1,2), sep='/')

master<-rbind(master, data.frame(dataID='ZAMO1',sp=p1$sp, colony=p1$Colony,
                                 trackID=p1$ID,date=trimws(p1$Date),
                                 time=trimws(p1$Time), latitude=p1$Latitude,longitude=p1$Longitude, breedstage=p1$breedstage))

# Nichol BRBO & RFBO

p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Nicholl_Chagos_BRBO.csv', h=T)

master<-rbind(master, data.frame(dataID='NICH1',sp='BRBO', colony='Danger',
                                 trackID=factor(p1$birdID),date=gsub('-', '/', substr(p1$date, 1,10)),
                                 time=substr(p1$date, 12, 19),
                                 latitude=p1$latitude,longitude=p1$longitude, breedstage='chick'))

p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Nicholl_Chagos_RFBO.csv', h=T)

p1$date<-ifelse(nchar(as.character(p1$date))<19, 
                paste0(as.character(p1$date), ':00'), as.character(p1$date))

master<-rbind(master, data.frame(dataID='NICH2',sp='RFBO', colony='Danger',
                                 trackID=factor(p1$birdID),date=gsub('-', '/', substr(p1$date, 1,10)),
                                 time=substr(p1$date, 12, 19),
                                 latitude=p1$latitude,longitude=p1$longitude, breedstage='chick'))

# Maxwell BRNO

p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Maxwell_Caribbean_BRNO.csv', h=T)

master<-rbind(master, data.frame(dataID='MAXW1',sp='BRNO', colony='Dry Tortugas',
                                 trackID=factor(p1$unique_id),
                                 date=paste('2016', sprintf("%02d", p1$month), sprintf("%02d", p1$day), sep='/'),
                                 time=paste(sprintf("%02d", p1$hour), sprintf("%02d", p1$min), sprintf("%02d", p1$sec), sep=':'),
                                 latitude=p1$lat,longitude=p1$lon,breedstage='chick'))

# Cecere WTSH

p1<-read.csv2('C:/seabirds/sourced_data/tracking_data/raw/Cecere_Seychelles_WTSH.csv', h=T)

master<-rbind(master, data.frame(dataID='CECE1',sp='WTSH', colony='Aride',
                                 trackID=factor(p1$RING),
                                 date=paste(substr(p1$DATA, 7,10),substr(p1$DATA, 4,5),substr(p1$DATA, 1,2), sep='/' ),
                                 time=paste(sprintf("%02d", p1$HOUR), sprintf("%02d", p1$MINUTE), sprintf("%02d", p1$SECOND), sep=':'),
                                 latitude=p1$N,longitude=p1$E, breedstage='chick'))

# Clarke BRBO, MABO, RFBO, LEFR, GRFR, RTTB 

p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Clarke_Ashmore_ALL.csv', h=T)
p1$colony='Adele'
p1[which(p1$CaptureSite=='Middle Island, Ashmore Reef'),]$colony<-'Mid Ashmore'
p1[which(p1$CaptureSite=="West Island, Ashmore Reef"),]$colony<-'West Ashmore'
p1$breedstage='Adele'
p1[which(p1$CaptureSite=='Middle Island, Ashmore Reef'),]$colony<-'Mid Ashmore'
p1$breedstage='chick'
p1[which(p1$BreedingStatus=='Incubating'),]$breedstage<-'incubation'
p1[which(p1$BreedingStatus=="EmptyNest"),]$breedstage<-'empty'

master<-rbind(master, data.frame(dataID='CLAR1',sp=p1$Species, colony=p1$colony,
                                 trackID=factor(p1$BandNo),
                                 date=gsub('-', '/', p1$Date),
                                 time=p1$Time,
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage=p1$breedstage))

# Mendez RFBO

setwd('C:/seabirds/sourced_data/tracking_data/raw/RFB-data-from-L.MENDEZ')

birds<-list.files(recursive=T)

counter=1
out<-NULL
for(i in 1:length(birds)){
 p1 <-read.csv2(birds[i], h=T)
 if(i>1){if(strsplit(birds[i], '/')[[1]][1] !=
            strsplit(birds[i-1], '/')[[1]][1]){counter=counter+1}}
 out<-rbind(out, 
       if("date.GMT._logger"%in% names(p1))
       {
         data.frame(dataID=paste0('MEND', counter),sp='RFBO', colony=strsplit(strsplit(birds[i], '/')[[1]][1], ' ')[[1]][1],
                    trackID=factor(p1$Fichier), date=p1$date.GMT._logger,time=p1$time.GMT._logger,
                    latitude=as.numeric(as.character(p1$lat)),longitude=as.numeric(as.character(p1$long)),
                    breedstage='chick')  
       }else{
          data.frame(dataID=paste0('MEND', counter),sp='RFBO', colony=strsplit(strsplit(birds[i], '/')[[1]][1], ' ')[[1]][1],
          trackID=if(!'Fichier'%in% names(p1)){strsplit(birds[i], '_')[[1]][1]}else{factor(p1$Fichier)},
          date=p1$Date,time=p1$Time,
          latitude=as.numeric(as.character(p1$Latitude)),longitude=as.numeric(as.character(p1$Longitude)),
          breedstage='chick')
          })
 print(i)
}

# sort space at the front of some times
out$time<-trimws(out$time)

# sort dates
out$date<-gsub('\\.', '/', out$date)

out$date<-paste(unlist(lapply(strsplit(out$date, '\\/'),
                    function(x){sprintf("%04d",as.numeric(x[3]))})),
      unlist(lapply(strsplit(out$date, '\\/'),
                            function(x){sprintf("%02d",as.numeric(x[2]))})),
      unlist(lapply(strsplit(out$date, '\\/'),
                    function(x){sprintf("%02d",as.numeric(x[1]))})), sep='/')
# fix to date '0014'
out[out$dataID=='MEND9',]$date<-paste0('2014', substr(out[out$dataID=='MEND9',]$date, 5, 10))
# fix reverse day/month
out[out$dataID=='MEND11',]$date<-paste(substr(out[out$dataID=='MEND11',]$date, 1, 4), substr(out[out$dataID=='MEND11',]$date, 9,10),
                                       substr(out[out$dataID=='MEND11',]$date, 6, 7), sep='/')
# 170 rows with blank lat long date time ID fields
out<-na.omit(out)
# sort ID for MEND12

out$trackID<-as.character(out$trackID)
out[out$dataID=='MEND12',]$trackID<-substr(out[out$dataID=='MEND12',]$trackID, 1, 6)

master<-rbind(master, out)

# Soanes data BRBO, MABO,SOTE

setwd('C:/seabirds/sourced_data/tracking_data/raw/Soanes')

soanes_sooty<-read.csv('Soanes_Caribbean_SOTN.csv')

master<-rbind(master, data.frame(dataID='SOAN1',sp='SOTE', colony='Dog',
                                 trackID=factor(paste('SOAN1', soanes_sooty$ID, sep='_')),
                                 date=paste(substr(soanes_sooty$Date,7,10), substr(soanes_sooty$Date,4,5), substr(soanes_sooty$Date,1,2), sep='/'),
                                 time=soanes_sooty$Time,
                                 latitude=soanes_sooty$y,longitude=soanes_sooty$x, breedstage='chick'))

# MABO all bird format
mabo1<-read_xlsx('Re__masked_booby_data_from_Anguilla/MB 2015.xlsx', sheet=1)
mabo1<-data.frame(mabo1, dataID='SOAN2');mabo1$ID=paste(mabo1$ID, mabo1$dataID, sep='_')
mabo2<-read.csv('Re__masked_booby_data_from_Anguilla/MB 2016 June .csv')
mabo2<-data.frame(mabo2, dataID='SOAN3');mabo2$ID=paste(mabo2$ID, mabo2$dataID, sep='_')
mabo2$Date<-paste(substr(mabo2$Date, 7,10), substr(mabo2$Date, 4,5), substr(mabo2$Date,1,2), sep='-')
mabo3<-read.csv('Re__masked_booby_data_from_Anguilla/MB april 2016.csv')
mabo3<-data.frame(mabo3, dataID='SOAN4');mabo3$ID=paste(mabo3$ID, mabo3$dataID, sep='_')
mabo3$Date<-paste(substr(mabo3$Date, 7,10), substr(mabo3$Date, 4,5), substr(mabo3$Date,1,2), sep='-')
mabo1<-rbind(mabo1[c(1:5, 7)], mabo2[,c(1:5, 7)], mabo3[,c(1:5, 7)])

mabo1$Time<-trimws(mabo1$Time)
mabo1$Date<-substr(mabo1$Date, 1, 10)

master<-rbind(master, data.frame(dataID=mabo1$dataID,sp='MABO', colony='Dog',
                                 trackID=factor(mabo1$ID),
                                 date=gsub('-', '/', mabo1$Date),
                                 time=mabo1$Time,latitude=mabo1$Latitude,
                                 longitude=mabo1$Longitude, breedstage='chick'))

# All Soanes data file by bird format
birds<-list.files(recursive=T)[-grep('zip',list.files(recursive=T))][1:276]
birds<-birds[-grep('gpx', birds)]


counter=5
soanes1<-NULL
for(i in 1:length(birds)){
  if(i>1){if(strsplit(birds[i], '/')[[1]][1] !=
             strsplit(birds[i-1], '/')[[1]][1]){counter=counter+1}}
  d1<-read.csv(birds[i])
  d1<-data.frame(dataID=paste0('SOAN', counter), sp='BRBO',
             colony='Dog',
             trackID=factor(substr(strsplit(birds[i], '/')[[1]][2], 1, nchar(strsplit(birds[i], '/')[[1]][2])-4)),
             date=d1$Date,time=substr(d1$Time,2,9),latitude=d1$Latitude,longitude=d1$Longitude,
             breedstage='chick') 
  soanes1<-rbind(soanes1, d1)}

soanes1$sp<-as.character(soanes1$sp)
soanes1$colony<-as.character(soanes1$colony)

soanes1[grep('MB', soanes1$trackID),]$sp<-'MABO'

soanes1[grep('Somb', soanes1$trackID),]$colony<-'Sombrero'
soanes1[grep('PPW', soanes1$trackID),]$colony<-'Prickly Pear'

master<-rbind(master, soanes1)

# Gilmour BRBO, MABO, RFBO, GRFR, MAFR
setwd('C:/seabirds/sourced_data/tracking_data/raw/Gilmour')

birds<-list.files(recursive=T)[-grep('zip',list.files(recursive=T))]
birds<-birds[grep('csv',birds)]
birds<-birds[63:length(birds)]#skip BFBO
birds<-birds[-grep('PB',birds)]# remove pena blanca from Diego G
birds<-birds[-grep('PJE',birds)]# remove pajereos from Diego G

counter=1
gilmour1<-NULL
for(i in 1:length(birds)){
  if(i>1){if(strsplit(birds[i], '/')[[1]][1] !=
             strsplit(birds[i-1], '/')[[1]][1]){counter=counter+1}}
  d1<-read.csv(birds[i])
  d1$Time<-as.character(d1$Time);d1$Date<-as.character(d1$Date)
  d1$Time<-trimws(d1$Time)
  d1$Date<-trimws(d1$Date)
  d1$Date<-ifelse(length(grep('-',d1$Date))>0, gsub('-', '/', d1$Date), d1$Date)
  d1$Date<-ifelse(nchar(as.character(d1$Date))<10, paste(strsplit(d1$Date[1], '\\/')[[1]][3],
                                sprintf("%02d", as.numeric(strsplit(d1$Date[1], '\\/')[[1]][1])),
                                sprintf("%02d", as.numeric(strsplit(d1$Date[1], '\\/')[[1]][2])), sep='/'), d1$Date)
  
  
  
  d1<-data.frame(dataID=paste0('GILM', counter), sp=substr(birds[i], 1, 4),
                 colony='Palmyra',
                 trackID=factor(substr(strsplit(birds[i], '/')[[1]][2], 1, nchar(strsplit(birds[i], '/')[[1]][2])-4)),
                 date=d1$Date,time=d1$Time,latitude=d1$Latitude,longitude=d1$Longitude,
                 breedstage='unknown') 
  gilmour1<-rbind(gilmour1, d1)}

gilmour1$breedstage<-as.character(gilmour1$breedstage)
gilmour1$colony<-as.character(gilmour1$colony)

gilmour1[grep('_inc_', gilmour1$trackID),]$breedstage<-'incubation'
gilmour1[grep('_brd_', gilmour1$trackID),]$breedstage<-'chick'

gilmour1[grep('_SJ_', gilmour1$trackID),]$colony<-'San Jorge'
gilmour1[grep('_CLR_', gilmour1$trackID),]$colony<-'Clarion'
gilmour1[grep('_TE_', gilmour1$trackID),]$colony<-'Tern'
gilmour1[grep('_AA_', gilmour1$trackID),]$colony<-'Alacranes'
gilmour1[grep('_II_', gilmour1$trackID),]$colony<-'Isabel'
gilmour1[grep('_IP_', gilmour1$trackID),]$colony<-'Pajaros'

table(gilmour1$sp, gilmour1$colony)

# Fix the same day issue!

gilmour1$date<-as.character(gilmour1$date)
gilmour1[which(is.na(gilmour1$date)),]$date<-'2014/07/01'

gilmour1$datetime <- as.POSIXct(strptime(paste(gilmour1$date, gilmour1$time, sep=""), "%Y/%m/%d %H:%M:%S"), "GMT")
gilmour1$tracktime <- as.double(gilmour1$datetime)

library(zoo)

gilmour1<-gilmour1%>%group_by(trackID)%>%
  mutate(dt=c(-999, diff(tracktime)))%>%as.data.frame()
gilmour1[gilmour1$dt> -1,]$dt<-NA
gilmour1$dt3<-na.locf(gilmour1$dt, na.rm = F)

gilmour2<-gilmour1%>%group_by(trackID)%>%
  mutate(dt4=seq(0, length(unique(dt3)), 1)[match(dt3, unique(dt3))])%>%
  mutate(date2=paste0(substr(date, 1,8), 
         sprintf("%02d", as.numeric(substr(date, 9,10))+dt4)))
# fixing impossible dates
gilmour2[substr(gilmour2$date2, 6,10) =='07/32',]$date2<-paste0(
  substr(gilmour2[substr(gilmour2$date2, 6,10) =='07/32',]$date2, 1, 5), '08/01')
gilmour2[substr(gilmour2$date2, 6,10) =='08/32',]$date2<-paste0(
  substr(gilmour2[substr(gilmour2$date2, 6,10) =='08/32',]$date2, 1, 5), '09/01')
gilmour2[substr(gilmour2$date2, 6,10) =='08/33',]$date2<-paste0(
  substr(gilmour2[substr(gilmour2$date2, 6,10) =='08/33',]$date2, 1, 5), '09/02')
gilmour2[substr(gilmour2$date2, 6,10) =='09/31',]$date2<-paste0(
  substr(gilmour2[substr(gilmour2$date2, 6,10) =='09/31',]$date2, 1, 5), '10/01')

gilmour2$datetime <- as.POSIXct(strptime(paste(gilmour2$date2, gilmour2$time, sep=""), "%Y/%m/%d %H:%M:%S"), "GMT")
gilmour2$tracktime <- as.double(gilmour2$datetime)

gilmour1$date<-gilmour2$date2
gilmour1$datetime<-NULL
gilmour1$tracktime<-NULL
gilmour1$d1<-NULL
gilmour1$dt3<-NULL
gilmour1$dt<-NULL

gilmour1<-gilmour1[gilmour1$trackID!='MABO_PAL_2017_inc_Female_1-1',] # rm crap trip

## back in

master<-rbind(master, gilmour1)

# Mcleay Crested Tern

p1<-NULL
for( i in 1:22)
{p1<-rbind(p1, data.frame(read_xlsx('C:/seabirds/sourced_data/tracking_data/raw/Mcleay_SouthAus_CRTE.xlsx', sheet=i,
                                   col_types=c('date', 'guess', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric')), 
                          ID=i))}

master<-rbind(master, data.frame(dataID='MCLE1',sp='CRTE', colony='Troubridge',
                                 trackID=factor(p1$ID),date=gsub('-', '/', as.character(p1$FixDate)), time=substr(p1$FixTime,1,8),
                                 latitude=p1$Lat,longitude=p1$Lon, breedstage='chick'))

# Dossa Senegalese Tern datasets

dos1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Dossa_Senegal_CATE.csv')

master<-rbind(master, data.frame(dataID='DOSA1',sp='CATE', colony='Senegal',
                                 trackID=factor(dos1$track_id),date=gsub('-', '/', as.character(dos1$date_gmt)),
                                 time=dos1$time_gmt,
                                 latitude=dos1$latitude,longitude=dos1$longitude, breedstage='incubation'))
dos1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Dossa_Senegal_ROTE1.csv')

master<-rbind(master, data.frame(dataID='DOSA2',sp='ROTE', colony='Senegal',
                                 trackID=factor(dos1$track_id),date=gsub('-', '/', as.character(dos1$date_gmt)),
                                 time=dos1$time_gmt,
                                 latitude=dos1$latitude,longitude=dos1$longitude, breedstage='incubation'))
dos1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/Dossa_Senegal_ROTE2.csv')

master<-rbind(master, data.frame(dataID='DOSA3',sp='ROTE', colony='Senegal',
                                 trackID=factor(dos1$track_id),date=gsub('-', '/', as.character(dos1$date_gmt)),
                                 time=dos1$time_gmt,
                                 latitude=dos1$latitude,longitude=dos1$longitude, breedstage='incubation'))

# Caroline Poli MABO
p1<-read.table('C:/seabirds/sourced_data/tracking_data/raw/Poli_Muertos_MABO.txt', sep='\t', h=T)

master<-rbind(master, data.frame(dataID='POLI1',sp='MABO', colony='Muertos',
                                 trackID=factor(p1$ID),date=gsub('-', '/', as.character(p1$DateGMT)),
                                 time=p1$TimeGMT,
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage='chick'))
# Julia Sommerfeld  MABO
p1<-read.table('C:/seabirds/sourced_data/tracking_data/raw/Sommerfeld_Phillip_MABO.txt', sep='\t', h=T)

master<-rbind(master, data.frame(dataID='SOMM1',sp='MABO', colony='Phillip',
                                 trackID=factor(p1$ID),date=paste(substr(p1$DateGMT, 7, 10), substr(p1$DateGMT, 4, 5),
                                                                  substr(p1$DateGMT, 1, 2), sep='/'),
                                 time=p1$TimeGMT,
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage='chick'))

# McDuie Swains MABO
lsf<-list.files('C:/seabirds/sourced_data/tracking_data/raw/Mcduie_Swins_MABO/SWAINS - 2013 deployments/SWAINS - 2013 deployments')
lsf<-lsf[grep('csv', lsf)]

p2<-NULL
for(i in 1:length(lsf))
{
  p1<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/raw/Mcduie_Swins_MABO/SWAINS - 2013 deployments/SWAINS - 2013 deployments/', lsf[i]))
  p1$trackID<-trimws(substr(lsf[i], 11,13))
  p2<-rbind(p2, p1)
}

master<-rbind(master, data.frame(dataID='MCDU1',sp='MABO', colony='Swains',
                                 trackID=factor(p2$trackID),date=p2$Date,
                                 time=p2$Time,
                                 latitude=p2$Latitude,longitude=p2$Longitude, breedstage='chick'))
# edit to change last track to BRBO
master[master$trackID=='y n' & !is.na(master$trackID),]$sp<-'BRBO'

# Austin RFBO, BRBO, MAFR

#brbo
p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/CaymanDataForMark/CaymanIslands_BrownBooby_GPSData_20162017(Austinetal).csv')
master<-rbind(master, data.frame(dataID='AUST1',sp='BRBO', colony='Cayman Brac',
                                 trackID=factor(p1$Bird_ID),date=paste(substr(p1$datetime_utc, 7, 10), substr(p1$datetime_utc, 4, 5),
                                substr(p1$datetime_utc, 1, 2), sep='/'),
                                 time=substr(p1$datetime_utc, 12, 19),
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage='chick'))
#rfbo GPS
p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/CaymanDataForMark/CaymanIslands_RedfootedBooby_GPSData_20162017(Austinetal).csv')
master<-rbind(master, data.frame(dataID='AUST2',sp='RFBO', colony='Little Cayman',
                                 trackID=factor(p1$Bird_ID),date=paste(substr(p1$datetime_utc, 7, 10), substr(p1$datetime_utc, 4, 5),
                                                                       substr(p1$datetime_utc, 1, 2), sep='/'),
                                 time=substr(p1$datetime_utc, 12, 19),
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage='chick'))
#rfbo gps-GSM
p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/CaymanDataForMark/CaymanIslands_RedfootedBooby_GSMData_20162017(Austinetal).csv')
master<-rbind(master, data.frame(dataID='AUST3',sp='RFBO', colony='Little Cayman',
                                 trackID=factor(p1$Bird_ID),date=paste(substr(p1$datetime_utc, 7, 10), substr(p1$datetime_utc, 4, 5),
                                                                       substr(p1$datetime_utc, 1, 2), sep='/'),
                                 time=substr(p1$datetime_utc, 12, 19),
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage='chick'))

#mag frigatebird
p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/CaymanDataForMark/CaymanIslands_Magnificent_Frigatebird_GPSData_2017.csv')

master<-rbind(master, data.frame(dataID='AUST4',sp='MAFR', colony='Little Cayman',
                                 trackID=factor(p1$Bird_ID),date=paste(substr(p1$datetime_utc, 7, 10), substr(p1$datetime_utc, 4, 5),
                                                                       substr(p1$datetime_utc, 1, 2), sep='/'),
                                 time=substr(p1$datetime_utc, 12, 19),
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage='chick'))
# St Helena tracking Oppel et al
p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/St Helena seabird tracking.csv')

#remove storm petrel
p1<-p1[p1$individual.taxon.canonical.name!='Oceanodroma castro',]
p1$individual.taxon.canonical.name<-as.character(p1$individual.taxon.canonical.name)
p1[p1$individual.taxon.canonical.name=='Phaethon aethereus',]$individual.taxon.canonical.name<-'RBTB'
p1[p1$individual.taxon.canonical.name=='Sula dactylatra',]$individual.taxon.canonical.name<-'MABO'
p1[p1$individual.taxon.canonical.name=='Anous stolidus',]$individual.taxon.canonical.name<-'BRNO'

p1$id<-'OPPE1'
p1[p1$individual.taxon.canonical.name=='BRNO',]$id<-'OPPE2'
p1[p1$individual.taxon.canonical.name=='RBTB',]$id<-'OPPE3'

master<-rbind(master, data.frame(dataID=p1$id,sp=p1$individual.taxon.canonical.name, colony='St Helena',
                                 trackID=factor(p1$tag.local.identifier),
                                 date=gsub('-', '/',substr(p1$timestamp, 1, 10)),
                                 time=substr(p1$timestamp, 12, 19),
                                 latitude=p1$location.lat,longitude=p1$location.long, breedstage='chick'))
# all breedstage as chick but need to update

# Gonzalez-Solis RBTB, BRBO

# RBTB Cabo verde multiple colonies
p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/GPS_PHAAET_J.Gonzalez-Solis_2017-2019_CaboVerde.csv')
p1$id<-'GONZ1'
p1[p1$Colony=='ICima',]$id<-'GONZ2'
p1[p1$Colony=='MLeao',]$id<-'GONZ3'
p1[p1$Colony=='PRoque',]$id<-'GONZ4'
p1[p1$Colony=='PSol_and_PRincon',]$id<-'GONZ5'
p1[p1$Colony=='Raso',]$id<-'GONZ6'
p1[p1$Colony=='SNegra',]$id<-'GONZ7'

p1$bs<-'chick'
p1[p1$Breed.Stage=='breeding',]$bs<-'unknown'
p1[p1$Breed.Stage=='incubation',]$bs<-'incubation'


master<-rbind(master, data.frame(dataID=p1$id,sp='RBTB', colony=p1$Colony,
                                 trackID=factor(p1$BirdId),
                                 date=paste(substr(p1$DateGMT, 7, 10), substr(p1$DateGMT, 4, 5),
                                            substr(p1$DateGMT, 1, 2), sep='/'),
                                 time=p1$TimeGMT,
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage=p1$bs))
# RBTB senegal

p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/GPS_PHAAET_J.Gonzalez-Solis_2015-2018_Senegal.csv')

p1$bs<-'chick'
p1[p1$Breed.Stage=='breeding',]$bs<-'unknown'
p1[p1$Breed.Stage=='incubating',]$bs<-'incubation'

master<-rbind(master, data.frame(dataID='GONZ8', sp='RBTB', colony=p1$Colony,
                                 trackID=factor(p1$BirdId),
                                 date=gsub('-', '/', p1$DateGMT),
                                 time=p1$TimeGMT,
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage=p1$bs))
# BRBO Cabo Verde

p1<-read.csv('C:/seabirds/sourced_data/tracking_data/raw/GPS_SULLEU_J.Gonzalez-Solis_2007-2018_CaboVerde.csv')

p1$bs<-'chick'
p1[p1$Breed.Stage=='breeding',]$bs<-'unknown'
p1[p1$Breed.Stage=='incubation',]$bs<-'incubation'

master<-rbind(master, data.frame(dataID='GONZ9', sp='BRBO', colony=p1$Colony,
                                 trackID=factor(p1$BirdId),
                                 date=gsub('-', '/', p1$DateGMT),
                                 time=p1$TimeGMT,
                                 latitude=p1$Latitude,longitude=p1$Longitude, breedstage=p1$bs))

# BLNO Cumming
a1<-st_read('C:/seabirds/sourced_data/tracking_data/raw/Cumming_BRNO/Noddy 063-68305 PP 2645 2019-12-13 10-38-34.kml')
a1$id='Nod68305'
a2<-st_read('C:/seabirds/sourced_data/tracking_data/raw/Cumming_BRNO/Noddy 063-68308 PP 2649 2019-12-13 21-22-39.kml')
a2$id='Nod68308'
a3<-st_read('C:/seabirds/sourced_data/tracking_data/raw/Cumming_BRNO/Noddy 063-68309 PP 2648 2019-12-14 23-15-31.kml')
a3$id='Nod68309'
a4<-st_read('C:/seabirds/sourced_data/tracking_data/raw/Cumming_BRNO/Noddy 063-68311 PP 2642 2019-12-13 10-42-01.kml')
a4$id='Nod68311'

p1<-rbind(a1,a2,a3,a4)

p1<-as.data.frame(p1)
p2<-strsplit(paste(p1$geometry), ',')
p2_out<-NULL
for(i in 1:length(p2))
{p2_out<-rbind(p2_out, data.frame(latitude=as.numeric(p2[[i]][2]),
                                  longitude=as.numeric(substr(p2[[i]][1], 3, nchar(p2[[i]][1])))))}

t1<-do.call(c, lapply(strsplit(as.character(p1$Name), ' '), function(x)x[2]))
t2<-do.call(c, lapply(strsplit(as.character(p1$Name), ' '), function(x)x[3]))

master<-rbind(master, data.frame(dataID='CUMM1',sp='BLNO', colony='Heron',
                                 trackID=p1$id, date=paste0('2019/12/',substr(t1, 1,2)),
                                 time=t2, p2_out, breedstage='chick'))


# temp write master, still some errors to fix
# remove all white from times
master$time<-trimws(master$time)
table(nchar(as.character(master$time)))

table(nchar(as.character(master$date)))
#master%>%group_by(dataID)%>%summarise_all(first)%>%print(n=60)

master$datetime <- as.POSIXct(strptime(paste(master$date, master$time, sep=""), "%Y/%m/%d %H:%M:%S"), "GMT")
master$tracktime <- as.double(master$datetime)

# rm duplicates
master<-master[which(!duplicated(paste0(master$trackID,master$tracktime))),]

master$datetime<-NULL

master<-master[-which(is.na(master$tracktime)),]

# make sure all tracks are temporally ordered 
master<-master[order(master$dataID, master$sp, master$colony,
                     master$trackID, master$tracktime, decreasing = F),]


write.csv(master, 'C:/seabirds/sourced_data/tracking_data/tracking_master.csv', quote=F, row.names=F)


# make summary table
 
master%>%group_by(sp)%>%summarise(nCol=length(unique(colony)),nBird=length(unique(trackID)),
                   nDset=length(unique(dataID)), own=paste(unique(dataID), collapse=','))


# make bbox for each dataset - note na.omit()
master_sf<-st_as_sf(na.omit(master), coords=c('longitude', 'latitude'), crs=4326)  #make sure the CRS is the same as your other data!

library(maps)
png('C:/seabirds/outputs/trackingmap.png', height=5, width=10, units='in', res=600 )

plot(latitude~longitude, data=master)
map('world', add=T, col=3)

dev.off()