# 12/12/19
# use GBR tracking data to generate candidate marine KBAs around 3 terrestrial KBAs 
# Trial using new mIBA (KBA) code from Oppel and Beal

#library(devtools)
#devtools::install_github("steffenoppel/track2iba", dependencies=TRUE) # add argument 'build_vignettes = FALSE' to speed it up
#not working
# use scripts seperately

library('sp')
library('sf')
library('smoothr')
library('raster')
library('tidyverse')
library('geosphere')
library('adehabitatHR')
library('foreach')
library('data.table')
library('purrr')
library('lubridate')
library('dplyr')
library('maptools')

source('C:/seabirds/sourced_data/code/new_mIBA/estSpaceUse.R')
source('C:/seabirds/sourced_data/code/new_mIBA/findKBA.R')
source('C:/seabirds/sourced_data/code/new_mIBA/findScale.R')
source('C:/seabirds/sourced_data/code/new_mIBA/repAssess.R')
source('C:/seabirds/sourced_data/code/new_mIBA/tripSummary.R')
source('C:/seabirds/sourced_data/code/new_mIBA/formatFields.R')
source('C:/seabirds/sourced_data/code/new_mIBA/tripSplit.R')


# read in data
wtsh<-read.csv('C:/seabirds/phd/analyses/tracking_data_pot/GPS_141516_clean_resamp_tripsplit_hmm.csv')
wtsh<-wtsh[wtsh$Colony=='Heron',]
wtsh<-wtsh[wtsh$Colony=='Heron',]

wtsh_meta<-read.csv('C:/seabirds/phd/analyses/tracking_data_pot/GPS_141516_trip_summary_attrib.csv')
wtsh<-wtsh[wtsh$trip_id %in% wtsh_meta[wtsh_meta$trip_type=='S',]$trip,]

brbo<-read.csv('C:/seabirds/phd/analyses/BRBO_raine/BRBO_raine_hmm.csv')

mabo1<-read.csv('C:/coral_fish/teaching/supervision/Charlotte_Dean/data/Swains_2014_tracking_trips_clean.csv')
mabo2<-read.csv('C:/coral_fish/teaching/supervision/Charlotte_Dean/data/Swains_2015_tracking_trips_clean.csv')

mabo2$Sex<-NULL; mabo2$Vel<-NULL
mabo<-rbind(mabo1,mabo2)
rm(mabo1);rm(mabo2)

# format data
brbo$ID<-factor(brbo$trip_id)
wtsh$ID<-factor(wtsh$trip_id)
mabo$ID<-factor(mabo$trip_id)

mabo<-mabo[mabo$ID!='-1',]

brbo<-formatFields(brbo, field_ID   = "ID", field_DateTime='DateTime',
                        field_Lon  = "Longitude", field_Lat  = "Latitude")
mabo<-formatFields(mabo, field_ID   = "ID", field_DateTime='DateTime',
                   field_Lon  = "Longitude", field_Lat  = "Latitude")
wtsh<-formatFields(wtsh, field_ID   = "ID", field_DateTime='DateTime',
                   field_Lon  = "Longitude", field_Lat  = "Latitude")
#colonies

raine=brbo %>% summarise(Longitude = first(Longitude), Latitude  = first(Latitude))
swains=mabo %>% summarise(Longitude = first(Longitude), Latitude  = first(Latitude))
heron=wtsh %>% summarise(Longitude = first(Longitude), Latitude  = first(Latitude))
                                                                           
                         
brbo_trips<-tripSplit(tracks = brbo, Colony= raine, 
  InnerBuff  = 3,ReturnBuff = 10, Duration   = 1,
  plotit     = TRUE,  rmColLocs  = F)

mabo_trips<-tripSplit(tracks = mabo, Colony= swains, 
                      InnerBuff  = 3,ReturnBuff = 10, Duration   = 1,
                      plotit     = TRUE,  rmColLocs  = F)

wtsh_trips<-tripSplit(tracks = wtsh, Colony= heron, 
                      InnerBuff  = 3,ReturnBuff = 10, Duration   = 1,
                      plotit     = TRUE,  rmColLocs  = F)


# brbo run

brbo_ts<-tripSummary(Trip=brbo_trips, Colony=raine)

Hvals<-findScale(brbo_trips, ARSscale = T, Trips_summary=brbo_ts, Colony=raine)

KDEs <- estSpaceUse( DataGroup = brbo_trips, Scale = 5, 
                    UDLev = 50, polyOut = TRUE, Res=2)

repr <- repAssess(brbo_trips, KDE = KDEs$KDE.Surface,
                  Scale = 5, Iteration = 100, BootTable = FALSE)

KBAs <- findKBA(KDE = KDEs, Represent = repr$out,
  UDLev = 50, Col.size = 2600, polyOut = TRUE,plotit = T) # col.size from batianoff

write_sf(KBAs, 'C:/seabirds/Birdlife_Aus/outputs/brbo_raine_kba.shp')

# mabo run

mabo_ts<-tripSummary(Trip=mabo_trips, Colony=swains)

Hvals<-findScale(mabo_trips, ARSscale = T, Trips_summary=mabo_ts, Colony=swains)

KDEs <- estSpaceUse( DataGroup = mabo_trips, Scale = 5, 
                     UDLev = 50, polyOut = TRUE, Res=2)

repr <- repAssess(mabo_trips, KDE = KDEs$KDE.Surface,
                  Scale = 5, Iteration = 1, BootTable = FALSE)

KBAs <- findKBA(KDE = KDEs, Represent = 80,
                UDLev = 50, Col.size = 100, polyOut = TRUE,plotit = T) # col.size from upper estimate Heatwole

write_sf(KBAs, 'C:/seabirds/Birdlife_Aus/outputs/mabo_swains_kba.shp')

# wtsh run

wtsh_ts<-tripSummary(Trip=wtsh_trips, Colony=heron)

Hvals<-findScale(wtsh_trips, ARSscale = T, Trips_summary=wtsh_ts, Colony=heron)

KDEs <- estSpaceUse( DataGroup = wtsh_trips, Scale = 15, 
                     UDLev = 50, polyOut = TRUE, Res=2)

repr <- repAssess(wtsh_trips, KDE = KDEs$KDE.Surface,
                  Scale = 15, Iteration = 100, BootTable = FALSE)

KBAs <- findKBA(KDE = KDEs, Represent = repr$out,
                UDLev = 50, Col.size = 26800, polyOut = TRUE,plotit = T) # col.size from Dyer

write_sf(KBAs, 'C:/seabirds/Birdlife_Aus/outputs/wtsh_heron_kba.shp')
