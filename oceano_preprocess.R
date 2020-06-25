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
r2<-raster(extent(r1), resolution=0.02,crs=CRS("+proj=longlat +datum=WGS84"), vals=1)
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

# read in embc attributed master

master_embc<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master_forage.csv')
# add colony column
master_embc$coly<-do.call(c, lapply(strsplit(as.character(master_embc$ID), '_'), function(x)x[3]))

# and trip quality control

t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')
t_qual<-t_qual[t_qual$manual_keep=='Y',]
# Fix yoda dates
t_qual$departure<-as.character(t_qual$departure)
t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$departure<-as.character(format(as.Date(as.POSIXlt(as.numeric(
  t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$departure),
  origin="1970-01-01", "GMT")),"%d/%m/%Y" ))

# and oceanographic month lookup

mo_look<-read.csv('C:/seabirds/data/dataID_month_lookup.csv')
# collapse by sp and colony
mo_look$coly<-do.call(c, lapply(strsplit(as.character(mo_look$ID), '_'), function(x)x[3]))
mo_look<-mo_look%>%filter(what=='decision')%>%group_by(sp, coly)%>%
summarise_at(vars(n1:n12), function(x){if('Y' %in% x){'Y'}else{''}})%>%as.data.frame()  

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


# 50% core area vs convex hull psuedo-abs approach

sp_groups<-
  list('BRBO', 'MABO', 'RFBO', 'WTSH', 'SOTE', 
       c('GRFR', 'LEFR', 'MAFR'), c('RBTB', 'RTTB'), 
       c('BRNO', 'LENO', 'BLNO'),c('CRTE', 'ROTE', 'CATE'))

names(sp_groups) <- c('BRBO', 'MABO', 'RFBO', 'WTSH', 'SOTE',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

hval<-0.03
  
for(k in sp_groups)
{
  dat<-master_embc[master_embc$sp %in% unlist(k),]
  
  # collapse colz within species
  colz_sp<-colz[colz$sp %in% unlist(k),]
  
  # Psuedo-absence sampled over all col dist accessibility within convex hull
  for(i in unique(dat$coly))
  {
  b1<-dat[dat$coly==i,]
  spdf<-SpatialPointsDataFrame(coords=b1[,c(7,8)],  data=data.frame(coly=b1$coly),
                               proj4string =CRS(projection(ex_templ)))
  r1<-crop(ex_templ, (extent(spdf)+0.3))
  r_pix<-as(r1,"SpatialPixels")
  #make dist
  
  r2<-rasterize(SpatialPoints(as(colz_sp[colz_sp$coly==i,],"Spatial")), r1)
  r2<-(distance(r2))/1000
  r2<-abs(r2-(max(values(r2)))) # leave vals in km
  
  hully<-st_as_sf(spdf)%>%summarise( geometry = st_combine( geometry ) ) %>%
    st_convex_hull()%>%st_buffer(dist=hval)
  #write_sf(st_as_sf(KDE.99), 'C:/seabirds/temp/brbo_for_50UD.shp')
  hully_ras<-rasterize(as(hully, 'Spatial'), r1) 
  hully_pts<-rasterToPoints(hully_ras, spatial=T)
  hully_pts$coly=i
  hully_pts$layer=0
  hully_pts$weight<-extract(r2, hully_pts)
  # Foraging within 50% UD
  bfor<-b1[b1$embc=='foraging',]
  
  spdf<-SpatialPointsDataFrame(coords=bfor[,c(7,8)],  data=data.frame(coly=bfor$coly),
                               proj4string =CRS(projection(ex_templ)))
  
  spdf$coly<-factor(spdf$coly)
  
  KDE.Surface <- kernelUD(spdf,same4all = F, h=hval, grid=r_pix)
  KDE.50 <- getverticeshr(KDE.Surface, percent = 50)
  #write_sf(st_as_sf(KDE.50), 'C:/seabirds/temp/brbo_for_50UD.shp')
  KDE.50_ras<-rasterize(KDE.50, r1) 
  KDE.50_pts<-rasterToPoints(KDE.50_ras, spatial=T)
  KDE.50_pts$coly=i
  KDE.50_pts$layer=1
  KDE.50_pts$weight=NA
  
  ###~~ Ocean data extract ~~###
  ext_pts<-rbind(st_as_sf(hully_pts),st_as_sf(KDE.50_pts))
  
  # lookup months
  moz<-which(mo_look[mo_look$sp %in% unlist(k) & mo_look$coly==i,3:14]=='Y')
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
 
  if(which(i==unique(dat$coly))==1){all_pts<-ext_pts}else{all_pts<-rbind(all_pts, ext_pts)}
  
    # bind up polygons for export
  names(KDE.50)[1]<-'coly'
  hully$coly<-i
  hully$PA=0
  KDE.50$area<-NULL
  KDE.50$PA=1
  stout2<-rbind(hully,st_as_sf(KDE.50 ))
  if(which(i==unique(dat$coly))==1){all_kerns<-stout2}else{all_kerns<-rbind(all_kerns, stout2)}
  
  print(i)
  
  }

#plot(all_pts[all_pts$coly=='Swains' & all_pts$layer==0, 'weight'])
# export polygons for gis
write_sf(all_kerns, paste0('C:/seabirds/data/GIS/', names(k), 'kernhull.shp'), delete_dsn=T)

write.csv(all_pts, paste0('C:/seabirds/data/modelling/kernhull_pts/', names(k), '_kernhull.csv'), quote=F, row.names=F)

print(k)

} #end multi-sp loop


# FAST kernel outputting to check errors

sp_groups<-
  list('BRBO', 'MABO', 'RFBO', 'WTSH', 'SOTE', 
       c('GRFR', 'LEFR', 'MAFR'), c('RBTB', 'RTTB'), 
       c('BRNO', 'LENO', 'BLNO'),c('CRTE', 'ROTE', 'CATE'))

names(sp_groups) <- c('BRBO', 'MABO', 'RFBO', 'WTSH', 'SOTE',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

hval<-0.03

k = sp_groups[1]

  dat<-master_embc[master_embc$sp %in% unlist(k),]
  # collapse colz within species
  colz_sp<-colz[colz$sp %in% unlist(k),]

  for(i in unique(dat$coly))
  {
    b1<-dat[dat$coly==i,]
    spdf<-SpatialPointsDataFrame(coords=b1[,c(7,8)],  data=data.frame(coly=b1$coly),
                                 proj4string =CRS(projection(ex_templ)))
    r1<-crop(ex_templ, (extent(spdf)+0.3))
    r_pix<-as(r1,"SpatialPixels")
 
    # Foraging within 50% UD
    bfor<-b1[b1$embc=='foraging',]
    spdf<-SpatialPointsDataFrame(coords=bfor[,c(7,8)],  data=data.frame(coly=bfor$coly),
                                 proj4string =CRS(projection(ex_templ)))
    
    spdf$coly<-factor(spdf$coly)
    
    KDE.Surface <- kernelUD(spdf,same4all = F, h=hval, grid=r_pix)
    KDE.50 <- getverticeshr(KDE.Surface, percent = 50)
    #write_sf(st_as_sf(KDE.50), 'C:/seabirds/temp/brbo_for_50UD.shp')

    # bind up polygons for export
    if(which(i==unique(dat$coly))==1){all_kerns<-st_as_sf(KDE.50)}else{all_kerns<-rbind(all_kerns, st_as_sf(KDE.50))}
    print(i)
  }

write_sf(all_kerns, paste0('C:/seabirds/temp/', names(k), 'kernhull_temp.shp'), delete_dsn=T)

##### OLD #####

# Create colony max range buffers and extract
# Pull in trip quality table
t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')
t_qual$auto_keep<-'Y'

# Booby trial
t_qual<-t_qual[grep('BRBO', t_qual$ID),]
t_qual[which(t_qual$duration>72),]$auto_keep<-'N'
t_qual$ID<-do.call(c, lapply(strsplit(as.character(t_qual$ID), '_'), function(x)x[3]))

# REMOVE DURATION == NA (NON RETURNS)

aggregate(max_dist~ID, t_qual, summary)
aggregate(max_dist~ID, t_qual[t_qual$auto_keep=='Y',], summary)

qplot(data=t_qual, x=max_dist, geom='histogram')+facet_wrap(~ID, scales='free')


# extract pred area
pred_a<-read_sf('C:/seabirds/data/GIS/pred_area.shp')
pred_a<-as(pred_a, 'Spatial')

# give pred_a pixel size of finest raster (bathy)
pred_ras<-rasterize(pred_a, bathy, field=1) # should attribute with hull rownumber if field in not specified
pred_pts<-rasterToPoints(pred_ras, spatial=T)

ex_sst<-extract(sst, pred_pts)
ex_chl<-extract(chl, pred_pts)
ex_front<-extract(fronts, pred_pts)
ex_bathy<-extract(bathy, pred_pts)
ex_slope<-extract(slope, pred_pts)

out3<-data.frame(pred_pts@coords, ex_sst, ex_chl, ex_front, ex_bathy, ex_slope)

out4<-na.omit(out3) # cut out land
write.csv(out4, 'C:/seabirds/data/pred_area_modelready.csv', quote=F, row.names=F)


# clip bathy to pred area extent as template for rasterize later on
#crop(bathy, pred_a, filename='C:/seabirds/data/GIS/pred_area_ras_template.tif')
