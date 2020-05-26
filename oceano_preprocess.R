# Oceanographic data pre-process and extract

# seasonal (4x .nc) CHL and SST import and create
# 1) annual climatology; 2) variabililty

library(raster)
library(ncdf4)

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

pfront_month<-raster('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_climatology_2006-2016.nc',
               varname='pfront')

pfront_sd<-calc(pfront_month, sd)

pfront_mn<-raster('C:/seabirds/sourced_data/pml_fronts/pml_CCI_SST_front-step3-sst_L3_tropics_1M_overall_monthly_2006-2016.nc',
               varname='pfront')


c_spring<-raster(chl_list[1],
                 varname='chlor_a')
c_summer<-raster('C:/seabirds/sourced_data/terra_chl_sst/T20001732019263.L3m_SCSU_CHL.x_chlor_a.nc',
                 varname='chlor_a')
c_summer<-raster('C:/seabirds/sourced_data/terra_chl_sst/T20002652019354.L3m_SCAU_CHL.x_chlor_a.nc',
                 varname='chlor_a')
c_summer<-raster('C:/seabirds/sourced_data/terra_chl_sst/T20003562019079.L3m_SCWI_CHL.x_chlor_a.nc',
                 varname='chlor_a')

st1<-stack(c_spring, c_summer, c_summer, c_summer)
c1<-calc(st1, sd)
c2<-calc(st1, mean)
