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

# read in embc attributed master

master_embc<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master_forage.csv')

# and trip quality control

t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')
t_qual<-t_qual[t_qual$manual_keep=='Y',]
# Fix yoda dates
t_qual$departure<-as.character(t_qual$departure)
t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$departure<-as.character(format(as.Date(as.POSIXlt(as.numeric(
  t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$departure),
  origin="1970-01-01", "GMT")),"%d/%m/%Y" ))

# attrib master

master_embc$sst<-NA
master_embc$sst_sd<-NA
master_embc$chl<-NA
master_embc$chl_sd<-NA
master_embc$mfr<-NA
master_embc$mfr_sd<-NA
master_embc$pfr<-NA
master_embc$pfr_sd<-NA
master_embc$bth<-NA
master_embc$slp<-NA

# we'll use the trip start month to extract data for that month for the
# whole trip, even if the trip spans multiple months.

# for dynamic vars
for (i in c('01','02','03','04','05','06', 
            '07','08','09','10', '11', '12'))
{
  # SST
  master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
                master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id, 9:10]<-
    
    extract(subset(sst_stack, c(as.numeric(i), 13)),
    master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
                master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id , 5:4])
  
  # CHL
  master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
            master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id, 11:12]<-
    
    extract(subset(chl_stack, c(as.numeric(i), 13)),
    master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
                master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id , 5:4])
  
  # mfront
  master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
                master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id, 13:14]<-
    
    extract(subset(mfront_stack, c(as.numeric(i), 13)),
    master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
                master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id , 5:4])
  
  # pfront
  master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
                master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id, 15:16]<-
    
    extract(subset(pfront_stack, c(as.numeric(i), 13)),
            master_embc[master_embc$ID %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$ID &
                          master_embc$trip_id %in% t_qual[which(substr(t_qual$departure, 4, 5)==i),]$trip_id , 5:4])
    
write.csv(master_embc, 'C:/seabirds/sourced_data/tracking_data/tracking_master_forage_extract.csv', 
          quote=F, row.names = F) 
print(i)
}

# for static vars

master_embc[,17]<-extract(bathy, master_embc[,5:4])
master_embc[,18]<-extract(slope, master_embc[,5:4])

write.csv(master_embc, 'C:/seabirds/sourced_data/tracking_data/tracking_master_forage_extract.csv', 
          quote=F, row.names = F)


# create and extract hulls

#create hulls around each dataID

hulls <- tracking %>%
  group_by( ID ) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_convex_hull()

#as sp:: objects for raster package

tracking<-as(tracking, 'Spatial')
hulls<-as(hulls, 'Spatial')
hulls$ID2<-1:length(hulls)

# give hulls pixels size of finest raster (bathy)
hull_ras<-rasterize(hulls, bathy) # should attribute with hull rownumber if field in not specified
hull_pts<-rasterToPoints(hull_ras, spatial=T)

ex_sst<-extract(sst, hull_pts)
ex_chl<-extract(chl, hull_pts)
ex_front<-extract(fronts, hull_pts)
ex_bathy<-extract(bathy, hull_pts)
ex_slope<-extract(slope, hull_pts)

hp<-left_join(as.data.frame(hull_pts@data), as.data.frame(hulls@data), by=c("layer"="ID2"))

out2<-data.frame(ID=hp$ID, ex_sst, ex_chl, ex_front, ex_bathy, ex_slope)

write.csv(out2, 'C:/seabirds/data/BRBO_hulls_modelready.csv', quote=F, row.names=F)

# 50% core area vs 99% UD psuedo-abs approach

brbo<-master_embc[master_embc$sp=='BRBO' & master_embc$trip_id %in% t_qual$trip_id,]
brbo$ID<-do.call(c, lapply(strsplit(as.character(brbo$ID), '_'), function(x)x[3]))
# same for colz
colz<-colz[grep('BRBO', colz$ID),]
colz$ID<-do.call(c, lapply(strsplit(as.character(colz$ID), '_'), function(x)x[3]))

## RM foraging points that are incorrectly identified as such e.g. sitting on island !

# Psuedo-absence sampled over all col dist accessibility within 99%UD


hval<-0.03

for(i in unique(brbo$ID))
{

b1<-brbo[brbo$ID==i,]
spdf<-SpatialPointsDataFrame(coords=b1[,c(5,4)],  data=data.frame(ID=b1$ID),
                             proj4string =CRS(projection(ex_templ)))
r1<-crop(ex_templ, (extent(spdf)+0.3))
r_pix<-as(r1,"SpatialPixels")
#make dist

r2<-rasterize(SpatialPoints(as(colz[colz$ID==i,],"Spatial")), r1)
r2<-(distance(r2))/1000
r2<-abs(r2-(max(values(r2)))) # leave vals in km

spdf$ID<-factor(spdf$ID)
KDE.Surface99 <- kernelUD(spdf,same4all = F, h=hval, grid=r_pix)
KDE.99 <- getverticeshr(KDE.Surface99, percent = 99 )
#write_sf(st_as_sf(KDE.99), 'C:/seabirds/temp/brbo_for_50UD.shp')
KDE.99_ras<-rasterize(KDE.99, r1) 
KDE.99_pts<-rasterToPoints(KDE.99_ras, spatial=T)
KDE.99_pts$ID=i
KDE.99_pts$layer=0
KDE.99_pts$weight<-extract(r2, KDE.99_pts)
# Foraging within 50% UD
bfor<-b1[b1$embc=='foraging',]

spdf<-SpatialPointsDataFrame(coords=bfor[,c(5,4)],  data=data.frame(ID=bfor$ID),
                             proj4string =CRS(projection(ex_templ)))

spdf$ID<-factor(spdf$ID)

KDE.Surface <- kernelUD(spdf,same4all = F, h=hval, grid=r_pix)
KDE.50 <- getverticeshr(KDE.Surface, percent = 50)
#write_sf(st_as_sf(KDE.50), 'C:/seabirds/temp/brbo_for_50UD.shp')
KDE.50_ras<-rasterize(KDE.50, r1) 
KDE.50_pts<-rasterToPoints(KDE.50_ras, spatial=T)
KDE.50_pts$ID=i
KDE.50_pts$layer=1
KDE.50_pts$weight=NA

stout<-rbind(st_as_sf(KDE.99_pts),st_as_sf(KDE.50_pts))
if(which(i==unique(brbo$ID))==1){all_pts<-stout}else{all_pts<-rbind(all_pts, stout)}

stout2<-rbind(st_as_sf(KDE.99 ),st_as_sf(KDE.50 ))
if(which(i==unique(brbo$ID))==1){all_kerns<-stout2}else{all_kerns<-rbind(all_kerns, stout2)}
print(i)
}

# vis and export to QGIS
plot(all_pts[all_pts$ID=='Swains' & all_pts$layer==0, 'weight'])

write_sf(all_kerns, 'C:/seabirds/temp/brbo_kerns_23Jun.shp')


ex_sst<-extract(sst, all_pts)
ex_chl<-extract(chl, all_pts)
ex_front<-extract(fronts, all_pts)
ex_bathy<-extract(bathy, all_pts)
ex_slope<-extract(slope, all_pts)

hp<-left_join(as.data.frame(KDE.50_pts@data), as.data.frame(KDE.50@data), by=c("layer"="ID2"))

out2<-data.frame(ID=hp$id, ex_sst, ex_chl, ex_front, ex_bathy, ex_slope)

write.csv(out2, 'C:/seabirds/data/BRBO_forKern_modelready.csv', quote=F, row.names=F)


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
