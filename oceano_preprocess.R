# Oceanographic data pre-process and extract

# seasonal (4x .nc) CHL and SST import and create
# 1) annual climatology; 2) variabililty

library(raster)
library(ncdf4)
library(dplyr)
library(sf)
library(adehabitatHR)

# Making annual mean and seasonal SD products of dynamic varibs and
# processing static varibs

# sst
sst_list<-list.files('C:/seabirds/sourced_data/terra_sst_month_clim', full.names=T)
sst_stack<-stack(sst_list, varname='sst4')

sst_mn<-calc(sst_stack, mean)
sst_sd<-calc(sst_stack, sd)

writeRaster(sst_mn, 'C:/seabirds/sourced_data/oceano_modelready/sst_mn.tif')
writeRaster(sst_sd, 'C:/seabirds/sourced_data/oceano_modelready/sst_sd.tif')

#chl
chl_list<-list.files('C:/seabirds/sourced_data/terra_chl_month_clim', full.names=T)
chl_stack<-stack(chl_list, varname='chl_ocx')

chl_stack<-log(chl_stack) # log prior to calc

chl_mn<-calc(chl_stack, mean)
chl_sd<-calc(chl_stack, sd)

writeRaster(chl_mn, 'C:/seabirds/sourced_data/oceano_modelready/chl_mn.tif')
writeRaster(chl_sd, 'C:/seabirds/sourced_data/oceano_modelready/chl_sd.tif')


#fronts
nc_open('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc')

# front strength
mfront_month<-stack('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc',
               varname='fronts_mean')

mfront_sd<-calc(mfront_month, sd)

mfront_mn<-raster('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_overall_monthly_2006-2016.nc',
               varname='fronts_mean')

writeRaster(mfront_mn, 'C:/seabirds/sourced_data/oceano_modelready/mfront_mn.tif')
writeRaster(mfront_sd, 'C:/seabirds/sourced_data/oceano_modelready/mfront_sd.tif')

# probability of front at pixel
pfront_month<-stack('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc',
                    varname='pfront')

pfront_sd<-calc(pfront_month, sd)

pfront_mn<-raster('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_overall_monthly_2006-2016.nc',
                  varname='pfront')

writeRaster(pfront_mn, 'C:/seabirds/sourced_data/oceano_modelready/pfront_mn.tif')
writeRaster(pfront_sd, 'C:/seabirds/sourced_data/oceano_modelready/pfront_sd.tif')

#bathy read and mosaic then make slope

b1<-raster('C:/seabirds/sourced_data/gebco/gebco_2020_n4.0_s-37.0_w-18.0_e173.0.nc')
b2<-raster('C:/seabirds/sourced_data/gebco/gebco_2020_n33.0_s-3.0_w-170.0_e-15.0.nc')
b3<-raster('C:/seabirds/sourced_data/gebco/gebco_2020_n27.0_s21.0_w120.0_e126.0.nc')

b4<-mosaic(b1, b2, b3, fun=min)
writeRaster(b4, 'C:/seabirds/sourced_data/oceano_modelready/bathy.tif')

# calc slope
terrain(b4, opt='slope', unit='degrees', neighbors=4, filename='C:/seabirds/sourced_data/oceano_modelready/slope.tif')

### Make standardised extraction grid @ 0.02 ~ 2km degree
r1<-raster('C:/seabirds/sourced_data/oceano_modelready/mfront_sd.tif')
r2<-raster(extent(r1), resolution=0.01,crs=CRS("+proj=longlat +datum=WGS84"), vals=1)
writeRaster(r2, 'C:/seabirds/sourced_data/oceano_modelready/extraction_template_1km.tif', overwrite=TRUE)


# read in EMbC classed tracking and extract oceano data

#read in monthly dynamic vars

# sst
sst_list<-list.files('C:/seabirds/sourced_data/terra_sst_month_clim', full.names=T)
sst_stack<-stack(sst_list, varname='sst4')
# layer 13 is sd product
sst_stack<-stack(sst_stack, 'C:/seabirds/sourced_data/oceano_modelready/sst_sd.tif')

#chl
chl_list<-list.files('C:/seabirds/sourced_data/terra_chl_month_clim', full.names=T)
chl_stack<-stack(chl_list, varname='chl_ocx')

chl_stack<-log(chl_stack) # log it
# layer 13 is sd product
chl_stack<-stack(chl_stack,'C:/seabirds/sourced_data/oceano_modelready/chl_sd.tif')

# front strength
mfront_stack<-stack('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc',
                    varname='fronts_mean')
# layer 13 is sd product
mfront_stack<-stack(mfront_stack, 'C:/seabirds/sourced_data/oceano_modelready/mfront_sd.tif')

# front probability
pfront_stack<-stack('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc',
                    varname='pfront')
# layer 13 is sd product
pfront_stack<-stack(pfront_stack, 'C:/seabirds/sourced_data/oceano_modelready/pfront_sd.tif')

bathy<-raster('C:/seabirds/sourced_data/oceano_modelready/bathy.tif')
slope<-raster('C:/seabirds/sourced_data/oceano_modelready/slope.tif')


# read in extraction template
ex_templ<-raster('C:/seabirds/sourced_data/oceano_modelready/extraction_template_1km.tif')

# read in col locs
colz<-st_read('C:/seabirds/data/GIS/trackingID_colony_locs.shp')
colz$sp<-do.call(c, lapply(strsplit(as.character(colz$ID), '_'), function(x)x[2]))
colz$coly<-do.call(c, lapply(strsplit(as.character(colz$ID), '_'), function(x)x[3]))
colz$spcol<-paste(colz$sp, colz$coly)

# load in trip quality control

t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')
t_qual<-t_qual[t_qual$manual_keep=='Y',]

# make WTSH l/S trip split
t_qual$WTSH_SL<-ifelse(t_qual$duration>(24*3) | t_qual$max_dist>300, 'L', 'S')
t_qual[grep('Aride', t_qual$ID),]$WTSH_SL<-'S'

# read in embc attributed master

master_embc<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master_forage.csv')
# add colony column
master_embc$coly<-do.call(c, lapply(strsplit(as.character(master_embc$ID), '_'), function(x)x[3]))
master_embc$spcol<-paste(master_embc$sp, master_embc$coly)

# EDIT WTSH species into short and long trip 'species' using t_qual
master_embc$sp<-as.character(master_embc$sp)
master_embc[master_embc$sp=='WTSH' & master_embc$trip_id %in% t_qual[t_qual$WTSH_SL=='S',]$trip_id,]$sp<-'WTST'
master_embc[master_embc$sp=='WTSH' & master_embc$trip_id %in% t_qual[t_qual$WTSH_SL=='L',]$trip_id,]$sp<-'WTLG'

# Secondary removal of trips or points that cause kernel errors (e.g. overnights on islands)
master_embc<-filter(master_embc, !(trip_id=='Nod683111' & embc=='foraging' & Latitude > '-23.2')) #overnight on island
master_embc<-filter(master_embc, !(spcol=='MABO Dog' & Latitude>18.18 & Latitude< 18.24 & Longitude > '-63.03' & Longitude < '-63.11')) # island deployment locs
master_embc<-filter(master_embc, !(spcol=='MABO Mid Ashmore' & Latitude> '-12.237' & Latitude< '-12.241' & Longitude > 122.977 & Longitude < 122.984)) # island nests
master_embc<-filter(master_embc, !(spcol=='RFBO Mid Ashmore' & Latitude> '-12.237' & Latitude< '-12.241' & Longitude > 122.977 & Longitude < 122.984))# island nests
master_embc<-filter(master_embc, !(spcol=='RFBO Mid Ashmore' & Latitude> '-12.259' & Latitude< '-12.262' & Longitude > 123.093 & Longitude < 123.099))# island nests
master_embc<-filter(master_embc, !(spcol=='LEFR Mid Ashmore' & Latitude> '-12.25' & Latitude< '-12.27' & Longitude > 123.02 & Longitude < 123.1))# island nests
master_embc<-filter(master_embc, !(spcol=='GRFR Mid Ashmore' & Latitude> '-12.240' & Latitude< '-12.247' & Longitude > 122.960 & Longitude < 123.972))# island nests
master_embc<-filter(master_embc, !(spcol=='BRBO Cayman Brac' & Latitude>19.67 & Latitude< 19.71 & Longitude > '-79.828' & Longitude < '-79.915'))# island nests
master_embc<-filter(master_embc, !(spcol=='BRBO Cayman Brac' & Latitude>19.69 & Latitude< 19.74 & Longitude > '-79.775' & Longitude < '-79.828'))# island nests
master_embc<-filter(master_embc, !(spcol=='BRBO Cayman Brac' & Latitude>19.71 & Latitude< 19.76 & Longitude > '-79.708' & Longitude < '-79.775'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Prickly Pear' & Latitude>18.17 & Latitude< 18.3 & Longitude > '-62.932' & Longitude < '-63.140'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Prickly Pear' & Latitude>18.17 & Latitude< 18.3 & Longitude > '-62.932' & Longitude < '-63.140'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Prickly Pear' & Latitude>18.284 & Latitude< 18.286 & Longitude > '-63.234' & Longitude < '-63.238'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Prickly Pear' & Latitude>17.61 & Latitude< 17.642 & Longitude > '-63.211' & Longitude < '-63.230'))# 
master_embc<-filter(master_embc, trip_id!='BB Somb 00940 M2') # rm deployment trip
master_embc<-filter(master_embc, !(spcol=='BRBO Dog' & Latitude>18.17 & Latitude< 18.3 & Longitude > '-62.932' & Longitude < '-63.140'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Dog' & Latitude>18.58 & Latitude< 18.6 & Longitude > '-63.42' & Longitude < '-63.432'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Dog' & Latitude>18.27 & Latitude< 18.28 & Longitude > '-63.273' & Longitude < '-63.278'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Dog' & Latitude>18.02 & Latitude< 18.1 & Longitude > '-63.077' & Longitude < '-63.167'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Pajarera' & Latitude>19.18 & Latitude< 19.25 & Longitude > '-104.542' & Longitude < '-104.722'))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Adele' & Latitude>'-19.91' & Latitude< '-19.94' & Longitude > 119.91 & Longitude < 119.925))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Adele' & Latitude>'-19.91' & Latitude< '-19.94' & Longitude > 119.91 & Longitude < 119.925))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Adele' & Latitude>'-16.84' & Latitude< '-16.87' & Longitude > 122.124 & Longitude < 122.154))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Adele' & Latitude>'-14.052' & Latitude< '-14.056' & Longitude > 121.775 & Longitude < 121.779))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Adele' & Latitude>'-15.057' & Latitude< '-15.068' & Longitude > 124.319 & Longitude < 124.334))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Adele' & Latitude>'-14.069' & Latitude< '-14.072' & Longitude > 125.771 & Longitude < 125.775))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Adele' & Latitude>'-12.672' & Latitude< '-12.673' & Longitude > 124.537 & Longitude < 124.540))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Mid Ashmore' & Latitude>'-12.237' & Latitude< '-12.240' & Longitude > 122.978 & Longitude < 122.983))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Mid Ashmore' & Latitude>'-10.980' & Latitude< '-10.984' & Longitude > 122.841 & Longitude < 122.845))# 
master_embc<-filter(master_embc, !(spcol=='BRBO Mid Ashmore' & Latitude>'-10.918' & Latitude< '-10.920' & Longitude > 123.028 & Longitude < 123.031))# 
#master_embc<-filter(master_embc, !(spcol=='BRBO Swains' & Latitude>'-21.8985' & Latitude< '-21.8991' & Longitude > 152.3694 & Longitude < 152.3701))# 
#master_embc<-filter(master_embc, !(spcol=='BRBO Swains' & Latitude>'-21.9797' & Latitude< '-21.9813' & Longitude > 152.4723 & Longitude < 152.4743))# 
#master_embc<-filter(master_embc, !(spcol=='BRBO Swains' & Latitude>'-21.9645' & Latitude< '-21.9663' & Longitude > 152.5678 & Longitude < 152.5687))# 


# load in oceanographic month lookup
mo_look<-read.csv('C:/seabirds/data/dataID_month_lookup.csv')

# collapse by sp and colony
mo_look$coly<-do.call(c, lapply(strsplit(as.character(mo_look$ID), '_'), function(x)x[3]))
mo_look$spcol<-paste(mo_look$sp, mo_look$coly)
mo_look<-mo_look%>%filter(what=='decision')%>%group_by(spcol)%>%
summarise_at(vars(n1:n12), function(x){if('Y' %in% x){'Y'}else{''}})%>%as.data.frame()  

# load in hval ref
hvals<-read.csv('C:/seabirds/data/dataID_hvals.csv') # same as updated code when using mag as reference scale
hvals_ref<-hvals%>%group_by(sp_group)%>%summarise(med_hval=median(mag, na.rm=T))
#p1<-ggplot(data=hvals, aes(x=sp_group, y=mag))+
#  geom_boxplot()+geom_point(aes(colour=factor(int_decision)), shape=1)+theme_bw()+
#  xlab('Species group')+ylab('log mean foraging range (km)')+scale_colour_discrete("GPS temporal \nresolution (s)")
#  png('C:/seabirds/outputs/kernel_scale.png',width = 6, height =6 , units ="in", res =600)
#  print(p1)
#  dev.off()


# export trips per month to decide extract time window for each dataset
mnth_lookup<-t_qual%>% filter(complete=='complete trip') %>%
          group_by(ID)%>%summarise(n1=length(which(substr(departure, 4,5)=='01')),
                                               n2=length(which(substr(departure, 4,5)=='02')),
                                               n3=length(which(substr(departure, 4,5)=='03')),
                                               n4=length(which(substr(departure, 4,5)=='04')),
                                               n5=length(which(substr(departure, 4,5)=='05')),
                                               n6=length(which(substr(departure, 4,5)=='06')),
                                               n7=length(which(substr(departure, 4,5)=='07')),
                                               n8=length(which(substr(departure, 4,5)=='08')),
                                               n9=length(which(substr(departure, 4,5)=='09')),
                                               n10=length(which(substr(departure, 4,5)=='10')),
                                               n11=length(which(substr(departure, 4,5)=='11')),
                                               n12=length(which(substr(departure, 4,5)=='12')))

d1<-as.data.frame(mnth_lookup)
d1[1:98, 2:13]<-''
d1.5<-as.data.frame(mnth_lookup)
d1$ID<-as.character(d1$ID)#
d1.5$ID<-as.character(d1.5$ID)
d1.5$what='month'
d1$what='decision'

d2<-rbind(d1.5,d1)
d2<-d2[order(d2$ID),]

#write.csv(d2, 'C:/seabirds/data/dataID_month_lookup.csv', quote=F, row.names=F) 

#### ~~~~ custom merging of WTSH Heron GPS (2015) with PTT kernels from 2011 and 2013 (Miller et al 2015, MEPS) ~~~~####
#kerns
gps_2015<-read_sf('C:/seabirds/data/GIS/WTLGkernhull.shp')%>%
  filter(spcol=='WTSH Heron' & PA==1)
ptt_2011<-read_sf('C:/seabirds/phd/analyses/paper2/spatial/LTHeronPTT2011.shp')%>%
  filter(id==25)%>% st_set_crs(st_crs(gps_2015)) #select 25% UD
ptt_2013<-read_sf('C:/seabirds/phd/analyses/paper2/spatial/LTHeronPTT2013.shp')%>%
  filter(id==25)%>% st_set_crs(st_crs(gps_2015))#select 25% UD
s1<-st_union(gps_2015,ptt_2013)
s2<-st_union(s1,ptt_2011)
s2<-s2%>%dplyr::select('spcol', 'PA')
#hull
gps_hull<-read_sf('C:/seabirds/data/GIS/WTLGkernhull.shp')%>%
  filter(spcol=='WTSH Heron' & PA==0)
ptt_pts<-read.csv('C:/seabirds/phd/analyses/paper2/spreads/heron_PTT_LT_only.csv')
ptt_pts<-ptt_pts[-grep('2012', ptt_pts$Date),]
ptt_pts<-ptt_pts[ptt_pts$Nest_id!=66,]
ptt_pts<-st_as_sf(ptt_pts,coords=c('Longitude', 'Latitude'), crs=4326)%>%
  summarise( geometry = st_combine( geometry ))
both_hulls<-st_union(gps_hull, ptt_pts)%>%st_convex_hull()

heron1<-rbind(s2, both_hulls)
wtlg_kerns<-filter(wtlg_kerns, spcol!='WTSH Heron')
wtlg_kerns<-rbind(wtlg_kerns, heron1)
write_sf(wtlg_kerns, 'C:/seabirds/data/GIS/WTLGkernhull.shp', delete_dsn=T)
####~~~ + ~~~~####

# 50% core area vs convex hull psuedo-abs approach

sp_groups<-
  list('BRBO', 'MABO', 'RFBO', 'SOTE', 'WTST', 'WTLG',
       c('GRFR', 'LEFR', 'MAFR'), c('RBTB', 'RTTB'), 
       c('BRNO', 'LENO', 'BLNO'),c('CRTE', 'ROTE', 'CATE'))

names(sp_groups) <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

for(m in 1:length(sp_groups))
{
  k<-sp_groups[m][[1]]
  k0<-k
  if(TRUE %in% (unlist(k)%in% c('WTST', 'WTLG'))){k0<-'WTSH'}
  
  myhv<-(hvals_ref[hvals_ref$sp_group==names(sp_groups[m]),]$med_hval)/111 # ~convert km to DD
  
  dat<-master_embc[master_embc$sp %in% unlist(k),]
  
  # collapse colz within species
  colz_sp<-colz[colz$sp %in% unlist(k0),]
  
  # Psuedo-absence sampled over all col dist accessibility within convex hull
  for(i in unique(dat$spcol))
  {
  b1<-dat[dat$spcol==i,]
  spdf<-SpatialPointsDataFrame(coords=b1[,c(7,8)],  data=data.frame(spcol=b1$spcol),
                               proj4string =CRS(projection(ex_templ)))
  r1<-crop(ex_templ, (extent(spdf)+0.3))
  r_pix<-as(r1,"SpatialPixels")
  #make dist
  
  r2<-rasterize(SpatialPoints(as(colz_sp[colz_sp$spcol==i,],"Spatial")), r1)
  r2<-(distance(r2))/1000
  r2<-abs(r2-(max(values(r2)))) # leave vals in km
  
  hully<-st_as_sf(spdf)%>%summarise( geometry = st_combine( geometry ) ) %>%
    st_convex_hull()%>%st_buffer(dist=myhv)
  #write_sf(st_as_sf(KDE.99), 'C:/seabirds/temp/brbo_for_50UD.shp')
  hully_ras<-rasterize(as(hully, 'Spatial'), r1) 
  hully_pts<-rasterToPoints(hully_ras, spatial=T)
  hully_pts$spcol=i
  hully_pts$layer=0
  hully_pts$weight<-extract(r2, hully_pts)
  # Foraging within 50% UD
  bfor<-b1[b1$embc=='foraging' & b1$ColDist>4000,] # remove 'halo' of foraging points at trip start/end
  if(substr(i, 1, 4)=='MAFR'|substr(i, 1, 4)=='GRFR'|substr(i, 1, 4)=='LEFR'|
     i=='CRTE Troubridge' | i=='CATE Senegal' | i=='ROTE Senegal' |
     i=='BRBO Palmyra'| i=='BRBO Swains'| i=='MABO Swains'){
    bfor<-b1[b1$embc=='relocating'& b1$ColDist>4000,]}
  
  spdf<-SpatialPointsDataFrame(coords=bfor[,c(7,8)],  data=data.frame(spcol=bfor$spcol),
                               proj4string =CRS(projection(ex_templ)))
  
  spdf$spcol<-factor(spdf$spcol)
  
  KDE.Surface <- kernelUD(spdf,same4all = F, h=myhv, grid=r_pix)
  KDE.50 <- getverticeshr(KDE.Surface, percent = 50)
  if(i=='BRBO Swains'){KDE.50 <-getverticeshr(KDE.Surface, percent = 65)} # more representative for small dataset
  #write_sf(st_as_sf(KDE.50), 'C:/seabirds/temp/brbo_for_50UD.shp')
  KDE.50_ras<-rasterize(KDE.50, r1) 
  KDE.50_pts<-rasterToPoints(KDE.50_ras, spatial=T)
  KDE.50_pts$spcol=i
  KDE.50_pts$layer=1
  KDE.50_pts$weight=NA
  
  ###~~ Ocean data extract ~~###
  ext_pts<-rbind(st_as_sf(hully_pts),st_as_sf(KDE.50_pts))
  
  # lookup months
  moz<-which(mo_look[mo_look$spcol==i,2:13]=='Y')
  # extract, 4 varibs use dynamic month lookup
  if(length(moz)>1){
    ext_pts$chl<-rowMeans(extract(subset(chl_stack, moz), ext_pts), na.rm=T)
    ext_pts$sst<-rowMeans(extract(subset(sst_stack, moz), ext_pts), na.rm=T)
    ext_pts$mfr<-rowMeans(extract(subset(mfront_stack, moz), ext_pts), na.rm=T)
    ext_pts$pfr<-rowMeans(extract(subset(pfront_stack, moz), ext_pts), na.rm=T)
    }else{
    ext_pts$chl<-extract(subset(chl_stack, moz), ext_pts)
    ext_pts$sst<-extract(subset(sst_stack, moz), ext_pts)
    ext_pts$mfr<-extract(subset(mfront_stack, moz), ext_pts)
    ext_pts$pfr<-extract(subset(pfront_stack, moz), ext_pts)}
  
  ext_pts$chl_sd<-extract(subset(chl_stack, 13), ext_pts)
  ext_pts$sst_sd<-extract(subset(sst_stack, 13), ext_pts)
  ext_pts$mfr_sd<-extract(subset(mfront_stack, 13), ext_pts)
  ext_pts$pfr_sd<-extract(subset(pfront_stack, 13), ext_pts)
  
  ext_pts$bth<-extract(bathy, ext_pts)
  ext_pts$slp<-extract(slope, ext_pts)
  
  # bind up for export as csv
  # export points to csv
  ext_pts$Longitude<-st_coordinates(ext_pts)[,1]
  ext_pts$Latitude<-st_coordinates(ext_pts)[,2]
  st_geometry(ext_pts)<-NULL
 
  if(which(i==unique(dat$spcol))==1){all_pts<-ext_pts}else{all_pts<-rbind(all_pts, ext_pts)}
  
    # bind up polygons for export
  names(KDE.50)[1]<-'spcol'
  hully$spcol<-i
  hully$PA=0
  KDE.50$area<-NULL
  KDE.50$PA=1
  stout2<-rbind(hully,st_as_sf(KDE.50 ))
  if(which(i==unique(dat$spcol))==1){all_kerns<-stout2}else{all_kerns<-rbind(all_kerns, stout2)}
  
  print(i)
  
  }

#plot(all_pts[all_pts$spcol=='Swains' & all_pts$layer==0, 'weight'])
# export polygons for gis
write_sf(all_kerns, paste0('C:/seabirds/data/GIS/', names(sp_groups[m]), 'kernhull.shp'), delete_dsn=T)

write.csv(all_pts, paste0('C:/seabirds/data/modelling/kernhull_pts/', names(sp_groups[m]), '_kernhull.csv'), quote=F, row.names=F)

print(k)

} #end multi-sp loop


# FAST kernel outputting to check errors

sp_groups<-
  list('BRBO', 'MABO', 'RFBO', 'SOTE', 'WTST', 'WTLG',
       c('GRFR', 'LEFR', 'MAFR'), c('RBTB', 'RTTB'), 
       c('BRNO', 'LENO', 'BLNO'),c('CRTE', 'ROTE', 'CATE'))

names(sp_groups) <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

for(m in 1:length(sp_groups))
{
  k<-sp_groups[m][[1]]
  k0<-k
  if(TRUE %in% (unlist(k)%in% c('WTST', 'WTLG'))){k0<-'WTSH'}
  
  myhv<-(hvals_ref[hvals_ref$sp_group==names(sp_groups[m]),]$med_hval)/111 # ~convert km to DD
  
  dat<-master_embc[master_embc$sp %in% unlist(k),]
  
  # collapse colz within species
  colz_sp<-colz[colz$sp %in% unlist(k0),]
  
  for(i in unique(dat$spcol))
  {
    b1<-dat[dat$spcol==i,]
   
    spdf<-SpatialPointsDataFrame(coords=b1[,c(7,8)],  data=data.frame(coly=b1$coly),
                                 proj4string =CRS(projection(ex_templ)))
    r1<-crop(ex_templ, (extent(spdf)+0.3))
    r_pix<-as(r1,"SpatialPixels")
 
    # Foraging within 50% UD
    bfor<-b1[b1$embc=='foraging' & b1$ColDist>4000,] # remove 'halo' of foraging points at trip start/end
    if(substr(i, 1, 4)=='MAFR'|
       i=='CRTE Troubridge' | i=='CATE Senegal' | i=='ROTE Senegal' |
       i=='BRBO Palmyra'| i=='BRBO Swains'| i=='MABO Swains'){bfor<-b1[b1$embc=='relocating'& b1$ColDist>4000,]}
    
    spdf<-SpatialPointsDataFrame(coords=bfor[,c(7,8)],  data=data.frame(coly=bfor$coly),
                                 proj4string =CRS(projection(ex_templ)))
    
    spdf$coly<-factor(spdf$coly)
    
    KDE.Surface <- kernelUD(spdf,same4all = F, h=myhv, grid=r_pix)
    KDE.50 <- getverticeshr(KDE.Surface, percent = 50)
    #rsuf <- as(KDE.Surface[[1]], "SpatialPointsDataFrame")
    #rsuf2<-rsuf[rsuf$ud>(max(rsuf$ud)*0.25),] #points have to be within ~ 0.75 UD
    #rsuf2<-rsuf2[sample(1:nrow(rsuf2), 
    #                    nrow(rsuf2[rsuf2$ud>(max(rsuf2$ud)*0.5),]),
    #                    replace=F, prob =rsuf2$ud),]
    #rsuf3<-rsuf[rsuf$ud>(max(rsuf$ud)*0.05),] #points have to be within ~ 0.95 UD
    #rsuf3<-rsuf3[sample(1:nrow(rsuf3), 
    #                    nrow(rsuf3[rsuf3$ud>(max(rsuf3$ud)*0.25),]),
    #                    replace=F, prob =rsuf3$ud),]
    #write_sf(st_as_sf(KDE.50), 'C:/seabirds/temp/brbo_for_50UD.shp')

    # bind up polygons for export
    if(which(i==unique(dat$spcol))==1){all_kerns<-st_as_sf(KDE.50)}else{all_kerns<-rbind(all_kerns, st_as_sf(KDE.50))}
    #if(which(i==unique(dat$spcol))==1){kerns75<-st_as_sf(rsuf2)}else{kerns75<-rbind(kerns75, st_as_sf(rsuf2))}
    #if(which(i==unique(dat$spcol))==1){kerns95<-st_as_sf(rsuf3)}else{kerns95<-rbind(kerns95, st_as_sf(rsuf3))}
    print(i)
  }

write_sf(all_kerns, paste0('C:/seabirds/temp/', names(sp_groups[m]), '50ud.shp'), delete_dsn=T)
#write_sf(kerns75, paste0('C:/seabirds/temp/', names(sp_groups[m]), 'kern75.shp'), delete_dsn=T)
#write_sf(kerns95, paste0('C:/seabirds/temp/', names(sp_groups[m]), 'kern95.shp'), delete_dsn=T)

}

##### Creation of prediction rasters #####

# One for GBR @ 2 km
# Another, globally @ 5 km

# extract pred area
tmpl2km<-raster('C:/seabirds/sourced_data/oceano_modelready/extraction_template_2km.tif')

pred_a<-read_sf('C:/seabirds/data/GIS/pred_area_large.shp')
pred_a<-as(pred_a, 'Spatial')

sst<-stack('C:/seabirds/sourced_data/oceano_modelready/sst_mn.tif',
           'C:/seabirds/sourced_data/oceano_modelready/sst_sd.tif')
chl<-stack('C:/seabirds/sourced_data/oceano_modelready/chl_mn.tif',
           'C:/seabirds/sourced_data/oceano_modelready/chl_sd.tif') # already logged
fronts<-stack('C:/seabirds/sourced_data/oceano_modelready/mfront_mn.tif',
              'C:/seabirds/sourced_data/oceano_modelready/mfront_sd.tif',
              'C:/seabirds/sourced_data/oceano_modelready/pfront_mn.tif',
              'C:/seabirds/sourced_data/oceano_modelready/pfront_sd.tif')
bathy<-raster('C:/seabirds/sourced_data/oceano_modelready/bathy.tif')
slope<-raster('C:/seabirds/sourced_data/oceano_modelready/slope.tif')
# global bathy and slope
#bathy<-raster('C:/seabirds/sourced_data/gebco_wgs84')
#slope<-terrain(bathy, opt='slope', unit='degrees', neighbors=4)


#stacks for wtsh Dec-Mar

# sst
sst_list<-list.files('C:/seabirds/sourced_data/terra_sst_month_clim', full.names=T)
sst_stack<-stack(sst_list[1:4], varname='sst4') 
#chl
chl_list<-list.files('C:/seabirds/sourced_data/terra_chl_month_clim', full.names=T)
chl_stack<-stack(chl_list[1:4], varname='chl_ocx')
chl_stack<-log(chl_stack) # log it

# front strength
mfront_stack<-stack('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc',
                    varname='fronts_mean')
mfront_stack<-subset(mfront_stack, 1:4)
# front probability
pfront_stack<-stack('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc',
                    varname='pfront')
pfront_stack<-subset(pfront_stack, 1:4)

# give pred_a pixel size of finest raster (2km)
pred_ras<-rasterize(pred_a, tmpl2km, field=1) # using 2 km rather than bathy
pred_pts<-rasterToPoints(pred_ras, spatial=T)

#pred_pts<-rasterToPoints(subset(fronts, 1), spatial=T) # make ~5km res global prediction layer

ex_sst<-extract(sst, pred_pts)
ex_chl<-extract(chl, pred_pts)
ex_front<-extract(fronts, pred_pts)
ex_bathy<-extract(bathy, pred_pts)
ex_slope<-extract(slope, pred_pts)

out3<-data.frame(pred_pts@coords, ex_sst, ex_chl, ex_front, ex_bathy, ex_slope)

out4<-na.omit(out3) # cut out land
write.csv(out4, 'C:/seabirds/data/pred_area_large_modelready_2km.csv', quote=F, row.names=F)
# write global 5km layer
#write.csv(out4, 'C:/seabirds/data/pred_area_global_modelready_5km.csv', quote=F, row.names=F)

#extract for wtsh only dyn_varibs
pt2<-SpatialPoints(out4[,c(1,2)], proj4string = CRS(projection(tmpl2km)))
ex_sst<-rowMeans(extract(sst_stack, pt2), na.rm=T)
ex_chl<-rowMeans(extract(chl_stack, pt2), na.rm=T)
ex_mfr<-rowMeans(extract(mfront_stack, pt2), na.rm=T)
ex_pfr<-rowMeans(extract(pfront_stack, pt2), na.rm=T)

out5<-data.frame(pt2@coords, ex_sst, ex_chl, ex_mfr, ex_pfr)

write.csv(out5, 'C:/seabirds/data/pred_area_large_modelready_2km_WTSHsummer.csv', quote=F, row.names=F)

# clip 2km pred_ras to pred area extent as template for rasterize later on
#crop(pred_ras, pred_a, filename='C:/seabirds/data/GIS/pred_area_large_ras_template.tif', overwrite=T)

#### Extract dynamic ocean varibs for each spcol within mean-maximum foraing rad @ 2km ####
# assumes env data is loaded and stacked as per kernel psuedo-abs extraction and month lookup is loaded

#degrade 1km ex_templ previously used to 2km resolution (save computation time)
r2<-raster(extent(ex_templ), resolution=0.02, crs=projection(ex_templ))

mo_look[mo_look$spcol=='SOTE chick',]$spcol<-'SOTE Rat'

# read in for global mean for rad data to get mean rad
for_rad<-read_sf('C:/seabirds/data/GIS/global_mean_foraging_radii.shp')
for_rad$spcol<-paste(for_rad$sp, for_rad$coly) # wedgies

# make dset for each sp-col combo inside radius
for( i in 1:nrow(for_rad))
{
my_rad<-for_rad[i,]
r1<-crop(r2, extent(my_rad))
rad_ras<-rasterize(as(my_rad, 'Spatial'), r1)
ext_pts<-st_as_sf(rasterToPoints(rad_ras, spatial=T))
###~~ Ocean data extract ~~###
# lookup months
spcol_lkup<-my_rad$spcol
if(substr(spcol_lkup, 1, 4)%in%c('WTST', 'WTLG')){
spcol_lkup<-paste('WTSH', substr(spcol_lkup,6, nchar(spcol_lkup)))}
moz<-which(mo_look[mo_look$spcol==spcol_lkup,2:13]=='Y')
# extract, 4 varibs use dynamic month lookup
if(length(moz)>1){
ext_pts$chl<-rowMeans(extract(subset(chl_stack, moz), ext_pts), na.rm=T)
ext_pts$sst<-rowMeans(extract(subset(sst_stack, moz), ext_pts), na.rm=T)
ext_pts$mfr<-rowMeans(extract(subset(mfront_stack, moz), ext_pts), na.rm=T)
ext_pts$pfr<-rowMeans(extract(subset(pfront_stack, moz), ext_pts), na.rm=T)
}else{
ext_pts$chl<-extract(subset(chl_stack, moz), ext_pts)
ext_pts$sst<-extract(subset(sst_stack, moz), ext_pts)
ext_pts$mfr<-extract(subset(mfront_stack, moz), ext_pts)
ext_pts$pfr<-extract(subset(pfront_stack, moz), ext_pts)}
ext_pts$chl_sd<-extract(subset(chl_stack, 13), ext_pts)
ext_pts$sst_sd<-extract(subset(sst_stack, 13), ext_pts)
ext_pts$mfr_sd<-extract(subset(mfront_stack, 13), ext_pts)
ext_pts$pfr_sd<-extract(subset(pfront_stack, 13), ext_pts)
ext_pts$bth<-extract(bathy, ext_pts)
ext_pts$slp<-extract(slope, ext_pts)
# bind up for export as csv
# export points to csv
ext_pts$Longitude<-st_coordinates(ext_pts)[,1]
ext_pts$Latitude<-st_coordinates(ext_pts)[,2]
st_geometry(ext_pts)<-NULL
ext_pts$layer<-my_rad$spcol
if(i==1){all_pts<-ext_pts}else{all_pts<-rbind(all_pts, ext_pts)}
print(i)
}

 
all_pts<-na.omit(all_pts)
all_pts[all_pts$bth>0,]$bth<-0 
all_pts$bth<-sqrt(all_pts$bth^2)# Remember nearshore front values which ==0 should be NA

names(all_pts)[12]<-'x'
names(all_pts)[13]<-'y'

write.csv(all_pts, 'C:/seabirds/data/global_col_mean_rad_env_2km_dynamic.csv', quote=F, row.names=F)
#write.csv(all_pts[,1:11], 'C:/seabirds/data/global_col_mean_rad_env_2km_dynamic_nocoord.csv', quote=F, row.names=F)
