# Oceanographic data pre-process and extract

# seasonal (4x .nc) CHL and SST import and create
# 1) annual climatology; 2) variabililty

library(raster)
library(ncdf4)
library(dplyr)
library(sf)


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

#read in
sst<-stack('C:/seabirds/sourced_data/oceano_modelready/sst_mn.tif',
      'C:/seabirds/sourced_data/oceano_modelready/sst_sd.tif')

chl<-stack('C:/seabirds/sourced_data/oceano_modelready/chl_mn.tif',
           'C:/seabirds/sourced_data/oceano_modelready/chl_sd.tif')

fronts<-stack('C:/seabirds/sourced_data/oceano_modelready/mfront_mn.tif',
             'C:/seabirds/sourced_data/oceano_modelready/mfront_sd.tif',
             'C:/seabirds/sourced_data/oceano_modelready/pfront_mn.tif',
             'C:/seabirds/sourced_data/oceano_modelready/pfront_sd.tif')

bathy<-raster('C:/seabirds/sourced_data/oceano_modelready/bathy.tif')
slope<-raster('C:/seabirds/sourced_data/oceano_modelready/slope.tif')

tracking<-read_sf('C:/seabirds/data/GIS/BRBO_foraging_25May.shp')

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