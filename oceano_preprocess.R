# Oceanographic data pre-process and extract

# seasonal (4x .nc) CHL and SST import and create
# 1) annual climatology; 2) variabililty

library(raster)
library(ncdf4)
library(dplyr)
library(sf)

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

# month to extract

t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')
t_qual<-t_qual[t_qual$manual_keep=='Y' & t_qual$complete=='complete trip',]

mnth_lookup<-t_qual%>%group_by(ID)%>%summarise(n1=length(which(substr(departure, 4,5)=='01')),
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

write.csv(mnth_lookup, 'C:/seabirds/data/dataID_month_lookup.csv')                                  


# compile tracking
l1<-list.files('C:/seabirds/sourced_data/tracking_data/foraging_embc')
for(i in l1)
{
  p1<-read.csv(paste0('C:/seabirds/sourced_data/tracking_data/foraging_embc/',
                      i))
  

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

# extract tracking
ex_sst<-extract(sst, tracking)
ex_chl<-extract(chl, tracking)
ex_front<-extract(fronts, tracking)
ex_bathy<-extract(bathy, tracking)
ex_slope<-extract(slope, tracking)

out<-data.frame(tracking@data, ex_sst, ex_chl, ex_front, ex_bathy, ex_slope)

write.csv(out, 'C:/seabirds/data/tracking_modelready.csv', quote=F, row.names=F)

# extract hulls
ex_sst<-extract(sst, hull_pts)
ex_chl<-extract(chl, hull_pts)
ex_front<-extract(fronts, hull_pts)
ex_bathy<-extract(bathy, hull_pts)
ex_slope<-extract(slope, hull_pts)

hp<-left_join(as.data.frame(hull_pts@data), as.data.frame(hulls@data), by=c("layer"="ID2"))

out2<-data.frame(ID=hp$ID, ex_sst, ex_chl, ex_front, ex_bathy, ex_slope)

write.csv(out2, 'C:/seabirds/data/BRBO_hulls_modelready.csv', quote=F, row.names=F)

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
