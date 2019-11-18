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