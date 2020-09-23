# Plotting and mapping
library(raster)
library(ggplot2)
library(ggspatial)
library(ggdendro)
library(dplyr)
#library(grid)
#library(gridExtra)
library(tidyr)
library(ggplotify)
library(patchwork)
library(vegan)
library(sf)
library(viridis)
library(geodist)
library(smoothr)


####~~~~ read in data ~~~~####
aucz_out<-read.csv('C:/seabirds/data/mod_validation_vals.csv')
matx_out<-read.csv('C:/seabirds/data/mod_clustering_vals.csv')
t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')
colz<-st_read('C:/seabirds/data/GIS/trackingID_colony_locs.shp')
spgroup_summ<-read.csv('C:/seabirds/data/sp_main_summary.csv')
gbr_cols<-read.csv('C:/seabirds/data/parks_gbr_colony_data.csv')

####~~~~ sp-col summary table ~~~~####

# formatting col
colz$sp<-do.call(c, lapply(strsplit(as.character(colz$ID), '_'), function(x)x[2]))
colz$coly<-do.call(c, lapply(strsplit(as.character(colz$ID), '_'), function(x)x[3]))
colz$spcol<-paste(colz$sp, colz$coly)
colz$Longitude<-st_coordinates(colz)[,1]
colz$Latitude<-st_coordinates(colz)[,2]
st_geometry(colz)<-NULL
colz<-colz%>%group_by(sp, coly)%>%summarise_all(first)

# formatting t_qual
t_qual<-t_qual[-which(t_qual$ID%in%c('MEND2_RFBO_Christmas', 'MEND3_RFBO_Christmas',
                                     'AUST3_RFBO_Little Cayman', 'GILM10_RFBO_Isabel')),]
t_qual$coly<-do.call(c, lapply(strsplit(as.character(t_qual$ID), '_'), function(x)x[3]))
t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$duration<-(as.numeric(as.character(t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$return))-
   as.numeric(as.character(t_qual[t_qual$ID=='YODA1_BRBO_Nakanokamishima',]$departure)))/3600 

#lookup col coords
t_qual<-left_join(t_qual, colz[,c(1,2,5,6)], by=c('sp', 'coly'))
# make WTSH l/S trip split
t_qual$WTSH_SL<-ifelse(t_qual$duration>(24*3) | t_qual$max_dist>300, 'L', 'S')
t_qual[grep('Aride', t_qual$ID),]$WTSH_SL<-'S'
t_qual$sp<-as.character(t_qual$sp)
t_qual[t_qual$sp=='WTSH' & t_qual$WTSH_SL=='S',]$sp<-'WTST'
t_qual[t_qual$sp=='WTSH' & t_qual$WTSH_SL=='L',]$sp<-'WTLG'
t_qual$sp_group<-t_qual$sp
t_qual[t_qual$sp %in% c('GRFR', 'LEFR', 'MAFR'),]$sp_group<-'FRBD'
t_qual[t_qual$sp %in% c('RBTB', 'RTTB'),]$sp_group<-'TRBD'
t_qual[t_qual$sp %in% c('BRNO', 'LENO', 'BLNO'),]$sp_group<-'NODD'
t_qual[t_qual$sp %in% c('CRTE', 'ROTE', 'CATE'),]$sp_group<-'TERN'

t_qual_unfilt<-t_qual

spcol_tab<-t_qual_unfilt%>%group_by(sp_group, sp, coly)%>%
  summarise(n_rawtracks=length(unique(trip_id)), col_long=first(Longitude), col_lat=first(Latitude))%>%as.data.frame()
# keep only gd trips
t_qual<-t_qual[t_qual$manual_keep=='Y',]
t_nused<-t_qual%>%group_by(sp_group, sp, coly)%>%
  summarise(n_usedtracks=length(unique(trip_id)))%>%as.data.frame()
spcol_tab<-data.frame(spcol_tab, n_usedtracks=t_nused$n_usedtracks)
# keep only returning trips
t_qual_ret<-t_qual[t_qual$complete=='complete trip',]
t_tripmetric<-t_qual_ret%>%filter(duration<296.1)%>%group_by(sp_group, sp, coly)%>%
  summarise(max_for=max(max_dist), med_for=median(max_dist), breedstage=paste(unique(breedstage), collapse=' '))%>%as.data.frame()
spcol_tab<-data.frame(spcol_tab, t_tripmetric[,c(4,5,6)])

# keep ALL trips as nobody using satellite GPS ie all loggers recovered from colonies
t_tripmetric2<-t_qual%>%filter(duration<296.1)%>%group_by(sp_group, sp, coly)%>%
  summarise(max_for_allret=max(max_dist), med_for_allret=median(max_dist))%>%as.data.frame()
spcol_tab<-data.frame(spcol_tab, t_tripmetric2[,c(4,5)])
# tidy
spcol_tab[spcol_tab$coly=='chick',]$coly<-'Rat' # edit to name
spcol_tab[spcol_tab$sp_group=='WTST' | spcol_tab$sp_group=='WTLG',]$sp<-'WTSH'
spcol_tab[spcol_tab$coly=='Adele' & spcol_tab$sp=='BRBO',]$max_for<-139.6  # manual edit
spcol_tab[spcol_tab$coly=='Adele' & spcol_tab$sp=='BRBO',]$max_for_allret<-139.6  # manual edit
spcol_tab[spcol_tab$coly=='Heron' & spcol_tab$sp_group=='WTLG',]$max_for_allret<-1150  # manual edit


#write.csv(spcol_tab, 'C:/seabirds/data/sp_col_summary.csv', quote=F, row.names=F)

# check max for ranges are not erroroneus
#tq_check<-t_qual_ret %>%filter(duration<296.1)%>% 
#     group_by(sp_group) %>%
#       filter(max_dist == max(max_dist))
#mst<-master_embc<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master_forage.csv')
#mst<-filter(mst, trip_id %in% tq_check$trip_id)
#write.csv(mst, 'C:/seabirds/temp/max_dist_trips.csv', quote=F, row.names=F)
#ggplot(data=mst, aes(x=Longitude, y=Latitude))+geom_point()+facet_wrap(~sp, scales='free')
# more detailed check

t_qual_ret_day<-t_qual[t_qual$complete=='complete trip',]
t_qual_ret_day<-t_qual_ret_day[t_qual_ret_day$breedstage %in% c('chick', 'incubation'),]
t_qual_ret_day$breedstage<-factor(t_qual_ret_day$breedstage)

t_qual_ret_day$day<-cut(t_qual_ret_day$duration, breaks=c(0,24*1:7, 999999),
                                   labels=(c('<1', '1-2', '2-3', '3-4','4-5', '5-6', '6-7', '>7')),right=F, include.lowest = T)
t_qual_ret_day$breedstage<-factor(t_qual_ret_day$breedstage)
t_tripmetric<-t_qual_ret_day%>%group_by(sp_group, sp, coly, breedstage, day)%>%
  summarise(ntrack=length(unique(trip_id)),max_for=max(max_dist))%>%as.data.frame()

#write.csv(t_tripmetric, 'C:/seabirds/data/sp_col_day_summary.csv', quote=F, row.names=F)
#### ~~~~ **** ~~~~ ####

####~~~~ sp-group summary table ~~~~####

spcol_tab<-read.csv('C:/seabirds/data/sp_col_summary.csv')

for_rang<-spcol_tab%>%group_by(sp_group)%>%summarise(min.for=min(max_for_allret),
                                                med.for=median(max_for_allret),
                                                max.for=max(max_for_allret))

sp_col_summr<-aucz_out%>%filter(as.character(spcol)!=as.character(Resample) & spcol!='MEAN')%>%
  group_by(sp, Resample)%>%
  summarise(mean_auc=mean(auc), sd_auc=sd(auc), max_auc=max(auc),
            mean_tss=mean(TSS), sd_tss=sd(TSS),max_tss=max(TSS))

sp_col_summr<-sp_col_summr%>%ungroup()%>%
  group_by(sp)%>%mutate(auc_rank=rank(-mean_auc, 'first'),tss_rank=rank(-mean_tss, 'first'))

mean_valz<-sp_col_summr%>%ungroup()%>%filter(!Resample %in% c('MultiCol', 'EnsembleRaw' , 'EnsembleNrm'))%>%
  group_by(sp)%>%
  summarise(mn_auc=mean(mean_auc), sd_auc=sd(mean_auc),
            mn_max_auc=mean(max_auc), sd_max_auc=sd(max_auc),
            mn_tss=mean(mean_tss), sd_tss=sd(mean_tss),
            mn_max_tss=mean(max_tss), sd_max_tss=sd(max_tss))

multi_rank<-sp_col_summr%>%filter(Resample=='MultiCol')%>%select(auc_rank, tss_rank)
ensemble_rank<-sp_col_summr%>%filter(Resample=='EnsembleRaw')%>%select(auc_rank, tss_rank)

#GBR pred
aucz_out%>%filter(as.character(spcol)!=as.character(Resample) & spcol%in%c('Raine', 'Swains', 'Heron')& auc>0.69)%>%group_by(sp, Resample)
GBR_pred<-sp_col_summr%>%ungroup()%>%filter(Resample!='MultiCol')%>%group_by(sp)

bind_out<-bind_cols(mean_valz, multi_rank, ensemble_rank,for_rang)

bind_out$mn_sumr<-paste0(round(bind_out$mn_auc, 2),'±', round(bind_out$sd_auc, 2))
bind_out$mx_sumr<-paste0(round(bind_out$mn_max_auc, 2),'±', round(bind_out$sd_max_auc, 2))
bind_out$mn_sumr_ts<-paste0(round(bind_out$mn_tss, 2),'±', round(bind_out$sd_tss, 2))
bind_out$mx_sumr_ts<-paste0(round(bind_out$mn_max_tss, 2),'±', round(bind_out$sd_max_tss, 2))


bind_out2<-bind_out%>%select(sp, mn_sumr, mx_sumr, mn_sumr_ts, mx_sumr_ts, auc_rank, tss_rank, auc_rank1, tss_rank1, min.for, med.for, max.for)
#write.csv(bind_out2, 'C:/seabirds/data/sp_main_summary.csv', quote=F, row.names=F)

#### ~~~~ **** ~~~~ ####

#### ~~~~ Make GBR colony radii min, med, max~~~~####
gbr_short<-gbr_cols%>%gather(species, trigger, -designation_name, -site_name, -designation_type, -Latitude, -Longitude )
gbr_short<-na.omit(gbr_short) # remove sp not at sites
gbr_short<-filter(gbr_short, !species %in% c('herald_petrel', 'silver_gull', 'australian_pelican',
                                   'bridled_tern', 'blacknaped_tern', 'little_tern',
                                   'roseate_tern', 'newcaledonianfairy_tern')) # rm unmodelled sp

#make lookup
gbr_short$mod_spgroup<-recode(gbr_short$species, brown_booby='BRBO',
                            masked_booby = 'MABO', redfooted_booby= 'RFBO',
                            lesser_frigatebird = 'FRBD', greater_frigatebird = 'FRBD',
                            redtailed_tropicbird = 'TRBD', wedgetailed_shearwater = 'WTST',
                            black_noddy='NODD', common_noddy='NODD',sooty_tern='SOTE',
                            caspian_tern='TERN', crested_tern='TERN', lessercrested_tern='TERN')

#write.csv(gbr_short%>%group_by(designation_name)%>%summarise(type=first(designation_type),
#                                                   site_names=paste(unique(site_name), collapse=' '),
#                                                   species=paste(unique(species), collapse=' '),
#                                                   spgroup=paste(unique(mod_spgroup), collapse=' ')),
#          'C:/seabirds/data/sites_bioregion_summary.csv', quote=F, row.names = F)

# filter and cols write as shapefile
# make spatial
#colz_sp<-gbr_short%>%st_as_sf(coords=c('Longitude', 'Latitude'), crs=4326)
#write_sf(colz_sp, 'C:/seabirds/data/GIS/parks_gbr_cols.shp')

# replicate data *3 for min, med and max buffers
gbr_rep<-bind_rows(gbr_short%>%mutate(rad_class='min'), 
                   gbr_short%>%mutate(rad_class='med'),
                   gbr_short%>%mutate(rad_class='max'))

#replicate wtsh long trips and bind back in
wtlong<-gbr_rep%>%filter(species=='wedgetailed_shearwater')
wtlong$mod_spgroup<-'WTLG'
gbr_rep<-bind_rows(gbr_rep, wtlong)%>%arrange(designation_name, site_name, mod_spgroup, species)

#do dist lookup
spgroup_dists<-spgroup_summ%>%select(Speces.group, Min.foraging.range, Median.foraging.range, Max.foraging.range)%>%
  gather(rad_class, dist, -Speces.group)
spgroup_dists$rad_class<-recode(spgroup_dists$rad_class,  Min.foraging.range='min',
                                Median.foraging.range ='med',  Max.foraging.range='max')
gbr_rep<-left_join(gbr_rep, spgroup_dists, by=c('mod_spgroup'='Speces.group', 'rad_class'))
# make spatial
gbr_rep<-gbr_rep%>%st_as_sf(coords=c('Longitude', 'Latitude'), crs=4326)

for(i in unique(gbr_rep$site_name))
{
  col<-gbr_rep%>%filter(site_name==i)
  colproj<-st_transform(col, crs=paste0('+proj=laea +lon_0=',st_coordinates(col)[1,1],
                                        '+lat_0=',st_coordinates(col)[2,2], '+ellps=WGS84'))
  colbuf<-colproj%>%st_buffer(dist=colproj$dist*1000)
  colbuf<-colbuf%>%st_transform(crs=4326)
  if(which(i==unique(gbr_rep$site_name))==1){buf_out<-colbuf}else{
  buf_out<-rbind(buf_out, colbuf)}
  print(i)
}

#write_sf(buf_out, 'C:/seabirds/data/GIS/foraging_radii.shp')

# not pred
# herald_petrel, silver_gull, australian_pelican, bridled_tern,
# blacknaped_tern, roseate_tern, newcaledonianfairy_tern

#### ~~~~ **** ~~~~ ####

#### ~~~~ Make GBR colony radii ALL distances ~~~~####
gbr_short<-gbr_cols%>%gather(species, trigger, -designation_name, -site_name, -designation_type, -Latitude, -Longitude )
gbr_short<-na.omit(gbr_short) # remove sp not at sites
gbr_short<-filter(gbr_short, !species %in% c('herald_petrel', 'silver_gull', 'australian_pelican',
                                             'bridled_tern', 'blacknaped_tern', 'little_tern',
                                             'roseate_tern', 'newcaledonianfairy_tern')) # rm unmodelled sp

#make lookup
gbr_short$mod_spgroup<-recode(gbr_short$species, brown_booby='BRBO',
                              masked_booby = 'MABO', redfooted_booby= 'RFBO',
                              lesser_frigatebird = 'FRBD', greater_frigatebird = 'FRBD',
                              redtailed_tropicbird = 'TRBD', wedgetailed_shearwater = 'WTST',
                              black_noddy='NODD', common_noddy='NODD',sooty_tern='SOTE',
                              caspian_tern='TERN', crested_tern='TERN', lessercrested_tern='TERN')


#replicate wtsh long trips and bind back in
wtlong<-gbr_short%>%filter(species=='wedgetailed_shearwater')
wtlong$mod_spgroup<-'WTLG'
gbr_short<-bind_rows(gbr_short, wtlong)

# gather species names within species groups
gbr_short<-gbr_short%>%group_by(designation_name, site_name, Latitude, Longitude, designation_type, mod_spgroup)%>%
  summarise(species=paste(species, collapse='_'), trigger=paste(trigger, collapse='_'))%>%
  arrange(designation_name, site_name, mod_spgroup, species)%>%as.data.frame()

#do dist lookup based on each radii
spcol_tab<-read.csv('C:/seabirds/data/sp_col_summary.csv')

gbr_rep<-left_join(gbr_short, spcol_tab[,c(1,9)], by=c('mod_spgroup'='Species.Group'))
# make spatial
gbr_rep<-gbr_rep%>%st_as_sf(coords=c('Longitude', 'Latitude'), crs=4326)

for(i in unique(gbr_rep$site_name))
{
  col<-gbr_rep%>%filter(site_name==i)
  colproj<-st_transform(col, crs=paste0('+proj=laea +lon_0=',st_coordinates(col)[1,1],
                                        '+lat_0=',st_coordinates(col)[2,2], '+ellps=WGS84'))
  colbuf<-colproj%>%st_buffer(dist=colproj$Max.roraging.range*1000)
  colbuf<-colbuf%>%st_transform(crs=4326)
  if(which(i==unique(gbr_rep$site_name))==1){buf_out<-colbuf}else{
    buf_out<-rbind(buf_out, colbuf)}
  print(i)
}

#write_sf(buf_out, 'C:/seabirds/data/GIS/foraging_radii_all_rad.shp')

#### ~~~~ **** ~~~~ ####

#### ~~~~ Make spgroup-colony foraging hotspots ~~~~ ####
all_rad<-read_sf('C:/seabirds/data/GIS/foraging_radii_all_rad.shp')

pred_list<-list.files('C:/seabirds/data/modelling/GBR_preds', full.names=T)
pred_list<-pred_list[grep('indivSUM.tif', pred_list)]
mod_pred<-stack(pred_list)

for( i in unique(all_rad$md_spgr))
{
  col<-all_rad[all_rad$md_spgr==i,]
  sp_ras1<-subset(mod_pred, paste0(i,'_indivSUM'))
  
  for(j in unique(col$site_nm))
  {
    col_sp<-filter(col, site_nm==j)%>%arrange(desc(Mx_rrg_))
    sp_ras<-crop(sp_ras1, extent(col_sp[1,]))
    
    # If long trip wtsh don't weight near colony as doesn't apply
    if(i=='WTLG'){sum_rad<-mask(sp_ras, as(col_sp[1,], 'Spatial'), updatevalue=0, updateNA=T)}else
      {
      for(k in 1:nrow(col_sp))
        {
          inrad<-mask(sp_ras, as(col_sp[k,], 'Spatial'), updatevalue=0, updateNA=T)
          ir2<-inrad
          ir2[ir2==0]<-NA
          q <- quantile(ir2, 0.7) # top xxx%
          toprad<-reclassify(inrad, c(-Inf, q, 0, q, Inf, 1))
          if(k==1){sum_rad<-toprad}else{sum_rad<-sum_rad+toprad}
          #plot(sum_rad)
        }
      }# close WTLG if loop
    # need to incorporate population information
    
    if(which(j==unique(col$site_nm))==1){mos_ras<-sum_rad}else{
      mos_ras <- mosaic(mos_ras, sum_rad, fun = max)} # mean works badly when overlap =0
    
    sr2<-sum_rad
    sr2[sr2==0]<-NA
    q2 <- quantile(sr2, 0.8) # top 20%
    hots<-reclassify(sum_rad, c(-Inf, q2, NA, q2, Inf, 1), right=F)
    s1<-st_as_sf(rasterToPolygons(hots, dissolve = T)) 
    s2<-drop_crumbs(s1, threshold=12000000)
    s3<-fill_holes(s2, threshold=54000000)
    #s4 <- smoothr::smooth(s3, method = "ksmooth", smoothness = 2)
    s3$dsgntn_n<-col_sp$dsgntn_n[1]
    s3$site_nm<-col_sp$site_nm[1]
    s3$dsgntn_t<-col_sp$dsgntn_t[1]
    s3$md_spgr<-col_sp$md_spgr[1]
    s3$species<-col_sp$species[1]
    s3$trigger<-col_sp$trigger[1]
    
    if(which(j==unique(col$site_nm))==1){core_pols<-s3}else{
      core_pols <- rbind(core_pols, s3)}
    
    plot(mos_ras)
    plot(core_pols, add=T)
    
  }
writeRaster(mos_ras, paste0('C:/seabirds/data/modelling/GBR_preds/col_radii_hotspots_', i, '.tif'), overwrite=T)
write_sf(core_pols, paste0('C:/seabirds/data/GIS/col_radii_core_hotspots_', i, '.shp'), delete_layer = T)
print(i)  
}

#### ~~~~ **** ~~~~ ####


#### ~~~~ Make foraging hotspot layer ~~~~ ####
pred_list<-list.files('C:/seabirds/data/modelling/GBR_preds', full.names=T)
pred_list<-pred_list[grep('indivSUM.tif', pred_list)]
mod_pred<-stack(pred_list)

for_rad<-read_sf('C:/seabirds/data/GIS/foraging_radii.shp')
rad_diss<-for_rad%>%group_by(md_spgr, rd_clss)%>%summarize(geometry = st_union(geometry))

for(i in unique(for_rad$md_spgr))
{
ext<-unlist(raster::extract(subset(mod_pred, paste0(i,'_indivSUM')),
                            as(filter(rad_diss, md_spgr==i & rd_clss=='max'), 'Spatial')))
mn<-min(ext, na.rm=T)
mx<-max(ext, na.rm=T)

r1<-subset(mod_pred, paste0(i,'_indivSUM'))
r1[values(r1)>mx]<-mx
r1[values(r1)<mn]<-mn

r2<-(r1-mn)/(mx-mn)# normalise (0-1)
#r3<-reclassify(r2, c(-Inf, 0.5, 0, 0.5,Inf,1))
if(i=='BRBO'){hotsp<-r2}else{hotsp<-hotsp+r2}
}
#writeRaster(hotsp, 'C:/seabirds/data/modelling/GBR_preds/hotspots.tif')

#### ~~~~ **** ~~~~ ####

#### ~~~~ Read in spatial data ~~~~ ####
for_rad<-read_sf('C:/seabirds/data/GIS/foraging_radii.shp')
# add sp colours
for_rad$md_spgr<-factor(for_rad$md_spgr, levels=c("BRBO", 'MABO', 'RFBO', 'FRBD', 'TRBD', 'WTST', "WTLG", "SOTE" , 'NODD', 'TERN'))
for_rad$spcol<-recode(for_rad$md_spgr, "BRBO"='#8dd3c7','MABO'='#ffffb3','RFBO'='#bebada','FRBD'='#fb8072', 'TRBD'='#80b1d3',
                      'WTST'='#fdb462','WTLG'='#b3de69',"SOTE"='#fccde5','NODD'='#d9d9d9','TERN'='#bc80bd')

land<-read_sf('C:/coral_fish/sourced_data/country_borders/TM_WORLD_BORDERS-0.3.shp')
gbr_reef<-read_sf('C:/seabirds/sourced_data/GBRMPA_Data Export/Great_Barrier_Reef_Features.shp')
gbrmp<-read_sf('C:/seabirds/sourced_data/GBRMPA_Data Export/GBRMP_BOUNDS.shp') 
colz<-read_sf('C:/seabirds/data/GIS/parks_gbr_cols.shp')
pred_list<-list.files('C:/seabirds/data/modelling/GBR_preds', full.names=T)
pred_list<-pred_list[grep('indivSUM.tif', pred_list)]
mod_pred<-stack(pred_list)
hotspots<-raster('C:/seabirds/data/modelling/GBR_preds/hotspots.tif')
# merge colz by md_spgr and rd_class for gbr-wide plots
rad_diss<-for_rad%>%group_by(md_spgr, rd_clss)%>%summarize(geometry = st_union(geometry))

#### ~~~~ GBR plot function ~~~~ ####
mk_gbrplot<-function(spg='TERN'){
  ext<-unlist(raster::extract(subset(mod_pred, paste0(spg,'_indivSUM')),
  as(filter(rad_diss, md_spgr==spg & rd_clss=='max'), 'Spatial')))
  mn<-min(ext, na.rm=T)
  mx<-max(ext, na.rm=T)
  
  r1<-subset(mod_pred, paste0(spg,'_indivSUM'))
  r1[values(r1)>mx]<-mx
  r1[values(r1)<mn]<-mn
  
  col_sp<-spg
  if(spg=='WTLG'){col_sp<-'WTST'} 

  p1<-ggplot() +
  layer_spatial(data=r1) +
  geom_sf(data=filter(colz, md_spgr==col_sp), shape = 23, fill = "darkred") +
  geom_sf(data=filter(rad_diss, md_spgr==spg), aes(colour=rd_clss), fill='NA') +
  geom_sf(data=gbrmp, col='white', fill='NA') +
  geom_sf(data=land, col='black', fill='grey') +
  theme_bw()+
  annotation_scale(location = "bl")+  
  annotation_north_arrow(location = "tr", which_north = "true")+
  labs(x='Longitude', y='Latitude')+
  coord_sf(xlim = c(142, 158), ylim = c(-29, -7), expand = FALSE)+
  scale_colour_manual('Forgaing radii', values=c('#00FFFF','#66FFCC', '#00FF66' ), labels=c(
    'Maximum', 'Median', 'Minimum'))+
  scale_fill_viridis('Likely\nforaging\nhabitat', limits=c(mn, mx), breaks=c(mn, mx), labels=c('low', 'high'),
                     option='magma', na.value = NA)

  return(p1)}

#theme(legend.position = c(0.2, 0.3), legend.background = element_rect(fill = "grey"),
 #     legend.key = element_rect(fill = "grey"))
# check higher res/different png dimensions
# patchwork <-p1 + p2 + p3 + p4 + plot_layout(ncol = 4,guides = "collect")

#### ~~~~ **** ~~~~ #####

#### ~~~~ Make GBR-wide plots ~~~~ ####
p_brbo<-mk_gbrplot(spg='BRBO')
p_mabo<-mk_gbrplot(spg='MABO')
p_rfbo<-mk_gbrplot(spg='RFBO')
p_frbd<-mk_gbrplot(spg='FRBD')
p_trbd<-mk_gbrplot(spg='TRBD')
p_wtst<-mk_gbrplot(spg='WTST')
p_wtlg<-mk_gbrplot(spg='WTLG')
p_sote<-mk_gbrplot(spg='SOTE')
p_nodd<-mk_gbrplot(spg='NODD')
p_tern<-mk_gbrplot(spg='TERN')

png(paste0('C:/seabirds/outputs/maps/gbr_wide/boobies2frigate.png'),width = 8.3, height =11.7 , units ="in", res =600)
p_brbo+ggtitle('A) Brown Booby')+p_mabo+ggtitle('B) Masked Booby')+
 p_rfbo+ggtitle('C) Red-footed Booby')+p_frbd+ggtitle('D) Frigatebird species-group')+
  plot_layout(ncol=2, nrow=2, guides = 'collect')&theme(legend.position = 'bottom')
dev.off()

png(paste0('C:/seabirds/outputs/maps/gbr_wide/tropbd2wtsh2sote.png'),width = 8.3, height =11.7 , units ="in", res =600)
p_trbd+ggtitle('E) Tropicbird species-group')+p_wtst+ggtitle('F) Wedge-tailed Shearwater short trips')+
    p_wtlg+ggtitle('G) Wedge-tailed Shearwater long trips')+p_sote+ggtitle('H) Sooty Tern')+
    plot_layout(ncol=2, nrow=2, guides = 'collect')&theme(legend.position = 'bottom')
dev.off()

png(paste0('C:/seabirds/outputs/maps/gbr_wide/nodd2tern.png'),width = 8.3, height =5.85 , units ="in", res =600)
p_nodd+ggtitle('I) Noddy species-group')+p_tern+ggtitle('J) Tern species-group')+
    plot_layout(ncol=2, guides = 'collect')&theme(legend.position = 'bottom')
dev.off()
#### ~~~~ **** ~~~~ #####

#### ~~~~ Make GBR-wide hotspot plot ~~~~ ####
mn<-min(values(hotspots), na.rm=T)
mx<-max(values(hotspots), na.rm=T)
hotspots<-(hotspots-mn)/(mx-mn)# normalise (0-1)
#hotspots_class<-cut(hotspots, breaks=seq(0, 1, 0.1))

p1<-ggplot() +
  layer_spatial(data=hotspots) +
  geom_sf(data=filter(colz), shape = 23, fill = "darkred") +
  geom_sf(data=gbrmp, col='white', fill='NA') +
  geom_sf(data=land, col='black', fill='grey') +
  theme_bw()+
  annotation_scale(location = "bl")+  
  annotation_north_arrow(location = "tr", which_north = "true")+
  labs(x='Longitude', y='Latitude')+
  coord_sf(xlim = c(142, 158), ylim = c(-29, -7), expand = FALSE)+
  scale_fill_viridis_b('Likely\nseabird\nforaging\nhabitat', option='magma',
                     breaks=c(seq(0.1, 0.9, 0.1)),labels=c('low', rep('', 7), 'high'), na.value = NA)

png(paste0('C:/seabirds/outputs/maps/gbr_wide/hotspots.png'),width = 8.3, height =5.85 , units ="in", res =600)
p1
dev.off()
#### ~~~~ **** ~~~~ #####

#### ~~~~ kba/site local plot function ~~~~ ####
mk_kbaplot<-function(site="Capricornia Cays KBA"){
  
  rad_diss_site<-filter(for_rad, dsgntn_n==site)%>%group_by(md_spgr, rd_clss)%>%
    summarize(spcol=first(spcol), species=first(species),geometry = st_union(geometry))
   plim<-st_bbox(filter(rad_diss_site,rd_clss=='med'))
  p1<-ggplot() +
    geom_sf(data=filter(gbr_reef, FEAT_NAME=='Reef'), fill='NA', colour=alpha('black',0.5))+
    geom_sf(data=filter(colz, dsgntn_n==site), shape = 23, fill = "darkred") +
    geom_sf(data=filter(rad_diss_site,rd_clss=='med'), aes(colour=spcol), fill='NA') +
    geom_sf(data=gbrmp, col='black', fill='NA') +
    geom_sf(data=land, col='black', fill='grey') +
    theme_bw()+
    annotation_scale(location = "bl")+  
    annotation_north_arrow(location = "tr", which_north = "true")+
    labs(x='Longitude', y='Latitude')+
    coord_sf(xlim = plim[c(1,3)], ylim = plim[c(2,4)], expand = T)+
    scale_colour_identity('Species',labels = unique(rad_diss_site$md_spgr),
        breaks = unique(rad_diss_site$spcol), guide = "legend")+ggtitle('A) Community')
  
  hotspot_local<-crop(hotspots, bbox(as(filter(rad_diss_site,rd_clss=='med'), 'Spatial')))
  p_last<-ggplot() +
    layer_spatial(data=hotspot_local)+
    geom_sf(data=filter(colz, dsgntn_n==site), shape = 23, fill = "darkred") +
    geom_sf(data=gbrmp, col='black', fill='NA') +
    geom_sf(data=land, col='black', fill='grey') +
    theme_bw()+
    annotation_scale(location = "bl")+  
    labs(x='Longitude', y='Latitude')+
    coord_sf(xlim = plim[c(1,3)], ylim = plim[c(2,4)], expand = F)+
      scale_fill_viridis(option='magma',na.value = NA)+
      ggtitle(paste0(LETTERS[length(unique(rad_diss_site$md_spgr))+2], ') Hotspots'))+
      theme(legend.position = 'none')
  
  p_spz<-list()
  for(i in unique(rad_diss_site$md_spgr))
  {
  r1<-subset(mod_pred, paste0(i,'_indivSUM'))
  beeb<-bbox(as(filter(rad_diss_site, md_spgr==i), 'Spatial'))
  beeb[c(1, 2)]<-beeb[c(1, 2)]-0.1
  beeb[c(3, 4)]<-beeb[c(3, 4)]+0.1
  r1<-crop(r1, beeb)
  mn<-min(values(r1), na.rm=T)+var(values(r1), na.rm=T) # small fudge
  mx<-max(values(r1), na.rm=T)-var(values(r1), na.rm=T)
  plim<-st_bbox(r1)  
  p2<-ggplot() +
    layer_spatial(data=r1) 
    if(i %in% c('BRBO' ,'MABO' ,'RFBO' ,'WTST', 'NODD','TERN')){
       p2<-p2+geom_sf(data=filter(gbr_reef, FEAT_NAME=='Reef'), fill='NA', colour=alpha('grey',0.5))}
    p2<-p2+geom_sf(data=filter(colz, dsgntn_n==site& md_spgr==i), shape = 23, fill = "darkred") +
    geom_sf(data=filter(rad_diss_site, md_spgr==i), aes(colour=rd_clss), fill='NA') +
    geom_sf(data=gbrmp, col='white', fill='NA') +
    geom_sf(data=land, col='black', fill='grey') +
    theme_bw()+
    annotation_scale(location = "bl")+  
    labs(x='Longitude', y='Latitude')+
    coord_sf(xlim = plim[c(1,3)], ylim = plim[c(2,4)], expand = F)+
    scale_colour_manual('Forgaing\nradii', values=c('#00FFFF','#66FFCC', '#00FF66' ), labels=c(
      'Max', 'Med', 'Min'))+
    scale_fill_viridis('Likely\nforaging\nhabitat',option='magma',
                       breaks=c(mn, mx), labels=c('low', 'high'),na.value = NA)
    
    if(1 %in% for_rad[for_rad$dsgntn_n==site & for_rad$md_spgr==i,]$trigger){
      p2<-p2+ggtitle(paste0(LETTERS[(which(i== unique(rad_diss_site$md_spgr)))+1], ') ', i, ' (Tr)'))
    } else{ 
    p2<-p2+ggtitle(paste0(LETTERS[(which(i== unique(rad_diss_site$md_spgr)))+1], ') ', i))}
    
  if(which(i== unique(rad_diss_site$md_spgr))!=1){p2<-p2+theme(legend.position = 'none')}
  p_spz[[which(i== unique(rad_diss_site$md_spgr))]]<-p2
  print(i)
  }
  
  pw_plot<-p1+p_spz+p_last+plot_layout(ncol=3, guides='collect')+
    plot_annotation(title = site)
  
  return(pw_plot)}
#### ~~~~ **** ~~~~ #####

#### ~~~~ Make local KBA  plots ~~~~ ####

for(k in unique(for_rad$dsgntn_n))
{
  multip<-mk_kbaplot(site=k) 
  png(paste0('C:/seabirds/outputs/maps/local_site/',gsub( ',', '',k),'.png'),width = 8.3, height =11.7 , units ="in", res =600)
  print(multip)
  dev.off()
  print(k)
}

#### ~~~~ **** ~~~~ #####

#### ~~~~ Make GBR validation plots ~~~~ ####
#
mkGBRval<-function(my.sp='BRBO', pred_col='Swains', with_reef=FALSE)
{
sp_val<-read.csv(paste0('C:/seabirds/data/modelling/GBR_validation/', my.sp, '_GBR_val.csv'))
# read foraging/hull shp
forhull<-read_sf(paste0('C:/seabirds/data/GIS/', my.sp, 'kernhull.shp'))
forhull$spcol<-substr(forhull$spcol, 6, nchar(forhull$spcol))
forhull<-filter(forhull, spcol==pred_col & PA==1)

if(my.sp=='BRBO'){gp1<-grep('MultiCol_',names(sp_val))
notin<-names(sp_val)[gp1][paste(my.sp, 'MultiCol', pred_col, sep='_')!=names(sp_val)[gp1]]
sp_val<-sp_val[,- which(notin==names(sp_val))]}

names(sp_val)[grep('MultiCol',names(sp_val))]<-paste0(my.sp, '_MultiCol')

sp_val<-filter(sp_val, spcol==pred_col)
#aucz_val<-filter(aucz_out, spcol==pred_col & sp==my.sp)
sp_val<-sp_val[,c(1, 2, 14:ncol(sp_val))]%>%tidyr::gather(model, pred, -spcol, -forbin, -Latitude, -Longitude)
sp_val$model<-substr(sp_val$model, 6, nchar(as.character(sp_val$model)))
sp_val$model<-gsub("\\.", " ", sp_val$model)
thresh_lkup<-sp_val%>%group_by(model)%>%
  summarise(thresh=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),
                          'best', best.method='youden', transpose=F)$threshold[1])
sp_val<-left_join(sp_val, thresh_lkup, by='model')
if('-Inf'%in%sp_val$thresh){sp_val[which(sp_val$thresh=='-Inf'),]$thresh<-0.5}
sp_val$PA<-ifelse(sp_val$pred>=sp_val$thresh, 1, 0)

sf_val<-sp_val%>%st_as_sf(coords=c('Longitude', 'Latitude'), crs=4326)
sf_val$model<-as.factor(sf_val$model)
sf_val$model<-factor(sf_val$model, levels=c(levels(sf_val$model)[levels(sf_val$model)!='MultiCol'], 'MultiCol'))
plim=st_bbox(sf_val)

p1<-ggplot()+geom_sf(data=sf_val, aes(colour=factor(PA)))

  if(with_reef==TRUE){p1<-p1+geom_sf(data=filter(gbr_reef, FEAT_NAME=='Reef'), fill='NA', colour=alpha('brown',0.2))}
  
p1<-p1+geom_sf(data=forhull, fill=NA)+facet_wrap(~model)+theme_bw()+
  coord_sf(xlim = plim[c(1,3)], ylim = plim[c(2,4)], expand = F)+
theme(legend.position = 'none')
return(p1)
}
#### ~~~~ **** ~~~~ #####

####~~~~ return GBR validation plots for GBR local sp ~~~~####
png('C:/seabirds/plots/brbo_swains_GBRvalidation2.png',width = 10, height =10 , units ="in", res =300)
print(mkGBRval(my.sp='BRBO', pred_col='Swains')) ;dev.off()

png('C:/seabirds/plots/brbo_raine_GBRvalidation2.png',width = 10, height =10 , units ="in", res =300)
print(mkGBRval(my.sp='BRBO', pred_col='Raine')) ;dev.off()

png('C:/seabirds/plots/mabo_swains_GBRvalidation2.png',width = 10, height =10 , units ="in", res =300)
print(mkGBRval(my.sp='MABO', pred_col='Swains')) ;dev.off()

png('C:/seabirds/plots/nodd_heron_GBRvalidation2.png',width = 10, height =10 , units ="in", res =300)
print(mkGBRval(my.sp='NODD', pred_col='Heron')) ;dev.off()

png('C:/seabirds/plots/wtst_heron_GBRvalidation2.png',width = 10, height =10 , units ="in", res =300)
print(mkGBRval(my.sp='WTST', pred_col='Heron')) ;dev.off()

png('C:/seabirds/plots/wtlg_heron_GBRvalidation2.png',width = 10, height =10 , units ="in", res =300)
print(mkGBRval(my.sp='WTLG', pred_col='Heron')) ;dev.off()


#### ~~~~ **** ~~~~ #####


####~~~~ AUC/TSS density plot ~~~~####

densdat<-aucz_out%>%filter(as.character(spcol)!=as.character(Resample) & spcol!='MEAN' &
                             !Resample %in% c('MultiCol', 'EnsembleRaw' , 'EnsembleNrm'))%>%
  select(sp, auc, TSS)%>%gather(variable, value, -sp)

densdat$sp<-factor(densdat$sp, levels=c('BRBO','MABO','RFBO','FRBD','TRBD','WTST','WTLG',
                                        'SOTE','NODD','TERN'))

p1<-ggplot(data=densdat, aes(x=value, stat(scaled), colour=variable))+
  geom_density()+facet_wrap(~sp)+theme_bw()+
  geom_vline(xintercept = 0.2, colour='coral2', linetype='dashed')+
  geom_vline(xintercept = 0.7, colour='seagreen3', linetype='dashed')+
  geom_hline(yintercept = 0)+
  scale_x_continuous(breaks=c(0,0.5, 1))+
  scale_colour_manual('Validation \nmetric',values=c( "seagreen3","coral2"), labels=c('AUC','TSS'))+
  theme(legend.position = c(0.8, 0.1),legend.justification = c(1, 0),
        panel.grid.minor.x = element_blank())+xlab('Metric value')+ylab('scaled density')

#png(paste0('C:/seabirds/outputs/sp_validation_comparison.png'),width = 6, height =6 , units ="in", res =600)
#p1
#dev.off()

#### ~~~~ **** ~~~~ ####

####~~~~ return cluster and heatmap validation plot function ~~~~####
mkVal<-function(my.sp='BRBO', my.metric='AUC', calc.niche=F)
{
  #https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html
  my.matx<-matx_out[matx_out$sp==my.sp,]
  my.aucz<-aucz_out[aucz_out$sp==my.sp,]
  
  if(my.metric=='AUC')
  {  
  d1<-with(my.matx[my.matx$Resample!=my.matx$spcol,c(8,9,3)], 
           structure(auc, Size = length(unique(my.matx$Resample)),
                     Labels = unique(my.matx$Resample),
                     Diag = F, Upper = FALSE,method = "user", class = "dist"))}
  if(my.metric=='TSS')
    {
    d1<-with(my.matx[my.matx$Resample!=my.matx$spcol,c(8,9,7)], 
             structure(TSS, Size = length(unique(my.matx$Resample)),
                       Labels = unique(my.matx$Resample),
                       Diag = F, Upper = FALSE,method = "user", class = "dist"))}
  
  d1<-(2-d1)
  hc1<-hclust(d1, method='average')
  hc_dend<-ggdendrogram(data = as.dendrogram(hc1), rotate = T)
  
  
  my.aucz$Resample<-as.character(my.aucz$Resample)
  my.aucz<-my.aucz[my.aucz$Resample!='EnsembleNrm',]
  my.aucz[my.aucz$Resample=='EnsembleRaw',]$Resample<-'Ensemble'
  my.aucz[my.aucz$Resample=='MultiCol',]$Resample<-'Multi-colony'
  my.aucz$Resample<-factor(my.aucz$Resample) 
  my.aucz$spcol<-factor(my.aucz$spcol)
  my.aucz$Resample<-factor(my.aucz$Resample, levels=c("Multi-colony", "Ensemble",paste(hc1$labels[hc1$order])))
  my.aucz$spcol<-factor(my.aucz$spcol,levels=c("MEAN", paste(hc1$labels[hc1$order])))
  
  if(my.metric=='AUC')
  {  
  val_plot<-ggplot(my.aucz, aes(x = spcol, y = Resample)) + 
    geom_raster(aes(fill=auc_bin)) +scale_fill_identity()+ 
    geom_text(aes(label=round(auc, 2)), size=2)+
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))+ylab('Predictions from')+xlab('Predicting to')}
  
  if(my.metric=='TSS')
  {
    val_plot<-ggplot(my.aucz, aes(x = spcol, y = Resample)) + 
      geom_raster(aes(fill=tss_bin)) +scale_fill_identity()+ 
      geom_text(aes(label=round(TSS, 2)), size=2)+
      theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11))+ylab('Predictions from')+xlab('Predicting to')}
  
  if(calc.niche==T){
  # niche overlap
  dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', my.sp, '_kernhull_sample.csv'))
  dat$X<-NULL
  enviro_std<-decostand(dat[,c(4:13)], method="standardize")
  enviro_rda<-rda(enviro_std, scale=T)
  #screeplot(enviro_rda)
  enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1))
  enviro.sites.scores<-data.frame(enviro.sites.scores,dat[,c(1,2)])
  enviro.species.scores<-as.data.frame(scores(enviro_rda, display='species'))
  enviro.species.scores$Predictors<-colnames(enviro_std)
  enviro.species.scores$spcol='X-Predictors'
  enviro.species.scores$PC1<-enviro.species.scores$PC1/30
  enviro.species.scores$PC2<-enviro.species.scores$PC2/30
  enviro.sites.scores$forbin<-relevel(enviro.sites.scores$forbin, ref='PsuedoA')
  
  pniche<-ggplot(data=enviro.sites.scores, aes(x=PC1, y=PC2))+
    geom_segment(data=enviro.species.scores, aes(y=0, x=0, yend=PC2, xend=PC1), arrow=arrow(length=unit(0.3,'lines')), colour='black')+
    geom_text(data=enviro.species.scores, aes(y=PC2, x=PC1, label=Predictors), colour='red', size=2)+
    geom_bin2d(aes(fill=forbin), alpha=0.6)+facet_wrap(~spcol)+theme_bw()+
    scale_fill_manual('Area',values=c( "seagreen3","coral2"), labels=c('Accessible \nhabitat','Core \nforaging'))
  
  return(list(hc1, val_plot, pniche))}else{
  
      return(list(hc1, val_plot))}
}
#### ~~~~ *** ~~~~ ####

#### ~~~~ AUC/TSS matrix and hclust plots ~~~~####
#https://cran.r-project.org/web/packages/ggplotify/vignettes/ggplotify.html

# BRBO 
brbo_auc<-mkVal('BRBO', 'AUC', calc.niche = F)
brbo_tss<-mkVal('BRBO', 'TSS', calc.niche = F)
brbo_auc_dend<-ggdendrogram(data = as.dendrogram(brbo_auc[[1]]), rotate = T)
brbo_tss_dend<-ggdendrogram(data = as.dendrogram(mabo_tss[[1]]), rotate = T)
# MABO
mabo_auc<-mkVal('MABO', 'AUC', calc.niche = F)
mabo_tss<-mkVal('MABO', 'TSS', calc.niche = F)
mabo_auc_dend<-ggdendrogram(data = as.dendrogram(brbo_auc[[1]]), rotate = T)
mabo_tss_dend<-ggdendrogram(data = as.dendrogram(mabo_tss[[1]]), rotate = T)
# RFBO
rfbo_auc<-mkVal('RFBO', 'AUC', calc.niche = F)
rfbo_tss<-mkVal('RFBO', 'TSS', calc.niche = F)
rfbo_auc_dend<-ggdendrogram(data = as.dendrogram(rfbo_auc[[1]]), rotate = T)
rfbo_tss_dend<-ggdendrogram(data = as.dendrogram(rfbo_tss[[1]]), rotate = T)
# FRBD
frbd_auc<-mkVal('FRBD', 'AUC', calc.niche = F)
frbd_tss<-mkVal('FRBD', 'TSS', calc.niche = F)
frbd_auc_dend<-ggdendrogram(data = as.dendrogram(frbd_auc[[1]]), rotate = T)
frbd_tss_dend<-ggdendrogram(data = as.dendrogram(frbd_tss[[1]]), rotate = T)
# TRBD
trbd_auc<-mkVal('TRBD', 'AUC', calc.niche = F)
trbd_tss<-mkVal('TRBD', 'TSS', calc.niche = F)
trbd_auc_dend<-ggdendrogram(data = as.dendrogram(trbd_auc[[1]]), rotate = T)
trbd_tss_dend<-ggdendrogram(data = as.dendrogram(trbd_tss[[1]]), rotate = T)
# WTST
wtst_auc<-mkVal('WTST', 'AUC', calc.niche = F)
wtst_tss<-mkVal('WTST', 'TSS', calc.niche = F)
wtst_auc_dend<-ggdendrogram(data = as.dendrogram(wtst_auc[[1]]), rotate = T)
wtst_tss_dend<-ggdendrogram(data = as.dendrogram(wtst_tss[[1]]), rotate = T)
# WTLG
wtlg_auc<-mkVal('WTLG', 'AUC', calc.niche = F)
wtlg_tss<-mkVal('WTLG', 'TSS', calc.niche = F)
wtlg_auc_dend<-ggdendrogram(data = as.dendrogram(wtlg_auc[[1]]), rotate = T)
wtlg_tss_dend<-ggdendrogram(data = as.dendrogram(wtlg_tss[[1]]), rotate = T)
# NODD
nodd_auc<-mkVal('NODD', 'AUC', calc.niche = )
nodd_tss<-mkVal('NODD', 'TSS', calc.niche = F)
nodd_auc_dend<-ggdendrogram(data = as.dendrogram(nodd_auc[[1]]), rotate = T)
nodd_tss_dend<-ggdendrogram(data = as.dendrogram(nodd_tss[[1]]), rotate = T)
# SOTE
sote_auc<-mkVal('SOTE', 'AUC', calc.niche = F)
sote_tss<-mkVal('SOTE', 'TSS', calc.niche = F)
sote_auc_dend<-ggdendrogram(data = as.dendrogram(sote_auc[[1]]), rotate = T)
sote_tss_dend<-ggdendrogram(data = as.dendrogram(sote_tss[[1]]), rotate = T)
# TERN
tern_auc<-mkVal('TERN', 'AUC', calc.niche = F)
tern_tss<-mkVal('TERN', 'TSS', calc.niche = F)
tern_auc_dend<-ggdendrogram(data = as.dendrogram(tern_auc[[1]]), rotate = T)
tern_tss_dend<-ggdendrogram(data = as.dendrogram(tern_tss[[1]]), rotate = T)

# AUC plots

brbomabo<-(brbo_auc[[2]]+ggtitle('A) Brown Booby'))+brbo_auc_dend+
  plot_layout(ncol=2, nrow=1, widths=c(3,1))
ggsave(plot=brbomabo, filename='C:/seabirds/data/modelling/plots/BRBO.png',width = 7, height =5)

maborfbo<-(rfbo_auc[[2]]+ggtitle('A) Masked Booby'))+mabo_auc_dend+(rfbo_auc[[2]]+ggtitle('B) Red-footed Booby'))+rfbo_auc_dend+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
ggsave(plot=maborfbo, filename='C:/seabirds/data/modelling/plots/MABO_RFBO.png',width = 7, height =10)

frbdtrbd<-(frbd_auc[[2]]+ggtitle('A) Frigatebirds'))+frbd_auc_dend+(trbd_auc[[2]]+ggtitle('B) Tropicbirds'))+trbd_auc_dend+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
ggsave(plot=frbdtrbd, filename='C:/seabirds/data/modelling/plots/FRBD_TRBD.png',width = 7, height =10)

wtstwtlg<-(wtst_auc[[2]]+ggtitle('A) Wedge-tailed Shearwater short trips'))+wtst_auc_dend+(wtlg_auc[[2]]+ggtitle('B) Wedge-tailed Shearwater long trips'))+wtlg_auc_dend+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
ggsave(plot=wtstwtlg, filename='C:/seabirds/data/modelling/plots/WTST_WTLG.png',width = 7, height =10)

noddsotetern<-(nodd_auc[[2]]+ggtitle('B) Noddies'))+nodd_auc_dend+
  (sote_auc[[2]]+ggtitle('B) Sooty Tern'))+sote_auc_dend+(tern_auc[[2]]+ggtitle('C) Terns'))+tern_auc_dend+
  plot_layout(ncol=2, nrow=3, widths=c(3,1))
ggsave(plot=noddsotetern, filename='C:/seabirds/data/modelling/plots/NODD_SOTE_TERN.png',width = 7, height =12)


# same for tss

brbomabo<-(brbo_tss[[2]]+ggtitle('A) Brown Booby'))+brbo_tss_dend+
  plot_layout(ncol=2, nrow=1, widths=c(3,1))
ggsave(plot=brbomabo, filename='C:/seabirds/data/modelling/plots/BRBO_tss.png',width = 7, height =5)

maborfbo<-(rfbo_tss[[2]]+ggtitle('A) Masked Booby'))+mabo_tss_dend+(rfbo_tss[[2]]+ggtitle('B) Red-footed Booby'))+rfbo_tss_dend+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
ggsave(plot=maborfbo, filename='C:/seabirds/data/modelling/plots/MABO_RFBO_tss.png',width = 7, height =10)

frbdtrbd<-(frbd_tss[[2]]+ggtitle('A) Frigatebirds'))+frbd_tss_dend+(trbd_tss[[2]]+ggtitle('B) Tropicbirds'))+trbd_tss_dend+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
ggsave(plot=frbdtrbd, filename='C:/seabirds/data/modelling/plots/FRBD_TRBD_tss.png',width = 7, height =10)

wtstwtlg<-(wtst_tss[[2]]+ggtitle('A) Wedge-tailed Shearwater short trips'))+wtst_tss_dend+(wtlg_tss[[2]]+ggtitle('B) Wedge-tailed Shearwater long trips'))+wtlg_tss_dend+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
ggsave(plot=wtstwtlg, filename='C:/seabirds/data/modelling/plots/WTST_WTLG_tss.png',width = 7, height =10)

noddsotetern<-(nodd_tss[[2]]+ggtitle('B) Noddies'))+nodd_tss_dend+
  (sote_tss[[2]]+ggtitle('B) Sooty Tern'))+sote_tss_dend+(tern_tss[[2]]+ggtitle('C) Terns'))+tern_tss_dend+
  plot_layout(ncol=2, nrow=3, widths=c(3,1))
ggsave(plot=noddsotetern, filename='C:/seabirds/data/modelling/plots/NODD_SOTE_TERN_tss.png',width = 7, height =12)



# If time
#library(rvg)
#library(officer)
#read_pptx('C:/coral_fish/outputs/portrait_template.pptx') %>%
#  add_slide(layout = "Title and Content", master = "Office Theme") %>%
#  ph_with(dml(ggobj=pw), location = ph_location_fullsize()) %>% 
#  print(target = 'C:/seabirds/temp/trial.pptx')





png(paste0('C:/seabirds/data/modelling/plots/temp.png'),width = 8, height =6 , units ="in", res =600)
grid.newpage()
print(brbo_auc[[2]], vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(brbo_auc[[1]], vp = viewport(x = 0.90, y = 0.58, width = 0.2, height = 0.82))
dev.off()
#### ~~~~ *** ~~~~ ####



####~~~~ local adaptation mantel and plots function ~~~~####
mklocada<-function(my.sp='BRBO', my.metric='AUC')
{
 
  #edit colony spatial info
  colz2<-colz
  colz2$sp<-do.call(c, lapply(strsplit(as.character(colz2$ID), '_'), function(x)x[2]))
  colz2$coly<-do.call(c, lapply(strsplit(as.character(colz2$ID), '_'), function(x)x[3]))
  colz2$spcol<-paste(colz2$sp, colz2$coly)
  colz2$Longitude<-st_coordinates(colz2)[,1]
  colz2$Latitude<-st_coordinates(colz2)[,2]
  st_geometry(colz2)<-NULL
  colz2<-colz2%>%group_by(sp, coly)%>%summarise_all(first)
  colz2[colz2$coly=='chick',]$coly<-'Rat'
  
  colz2$sp_group<-colz2$sp
  colz2[colz2$sp %in% c('GRFR', 'LEFR', 'MAFR'),]$sp_group<-'FRBD'
  colz2[colz2$sp %in% c('RBTB', 'RTTB'),]$sp_group<-'TRBD'
  colz2[colz2$sp %in% c('BRNO', 'LENO', 'BLNO'),]$sp_group<-'NODD'
  colz2[colz2$sp %in% c('CRTE', 'ROTE', 'CATE'),]$sp_group<-'TERN'
  colz2
  
  my.matx<-matx_out[matx_out$sp==my.sp,]
  
  if(my.metric=='AUC')
  {  
    d_pred<-with(my.matx[my.matx$Resample!=my.matx$spcol,c(10,11,5)], 
             structure(auc, Size = length(unique(my.matx$Resample)),
                       Labels = unique(my.matx$Resample),
                       Diag = F, Upper = FALSE,method = "user", class = "dist"))}
  if(my.metric=='TSS')
  {
    d_pred<-with(my.matx[my.matx$Resample!=my.matx$spcol,c(10,11,9)], 
             structure(TSS, Size = length(unique(my.matx$Resample)),
                       Labels = unique(my.matx$Resample),
                       Diag = F, Upper = FALSE,method = "user", class = "dist"))}
    
    # geographic dist
  
  col_lkup<-unique(my.matx$spcol)
  if(my.sp=='TERN' |my.sp=='FRBD'){col_lkup<-substr(col_lkup, 6, nchar(as.character(col_lkup)))}
  
    gd1<-filter(colz2, coly %in% col_lkup)%>%group_by(sp_group, sp, coly)%>%
      summarise_all(first)%>%as.data.frame()
    
    my.sp2<-my.sp
    if(my.sp=='WTST'|my.sp=='WTLG'){my.sp2<-'WTSH'}
    gd1<-filter(gd1, sp_group==my.sp2)
    
    #row.names(gd1)<-gd1$coly
    gd1<-gd1[c('Longitude', 'Latitude')]
    
    geo1<-geodist(gd1)
    df <- data.frame( Pred.perf=as.vector(d_pred), 
                      Geo.dist=geo1[lower.tri(geo1)])
    
    mn<-mantel(d_pred, as.dist(geo1), permutations=9999)
    
    p1<-ggplot(df,aes(x=Geo.dist,y=Pred.perf)) + geom_point(shape=1)+
      xlab('Pairwise colony geographic distance')+ylab('Pairwise colony predictive performance')+
      annotate(geom = 'text', label = paste('Mantel r =', round(mn$statistic, 3), 'p =', round(mn$signif, 3)),
               x = Inf, y = Inf, hjust = 1, vjust = 1)+theme_bw()
 
    # env niche
    dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', my.sp, '_kernhull_sample.csv'))
    dat$X<-NULL
    if(my.sp=='RFBO'){dat<-dat[dat$spcol!='Christmas',]}
    med_env<-dat[,c(2,4:13)]%>%group_by(spcol)%>%summarise_all(median, na.rm=T)%>%as.data.frame
    row.names(med_env)<-med_env$spcol
    env_std<-decostand(med_env[,c(2:11)], method="standardize")
    d_env<-dist(env_std, 'euclidean')
    
    df <- data.frame( Pred.perf=as.vector(d_pred), 
                      Ocean.cond=as.vector(d_env))
    
    mn<-mantel(d_pred, d_env, permutations=9999)
    
    p2<-ggplot(df,aes(x=Ocean.cond,y=Pred.perf)) + geom_point(shape=1)+
      xlab('Pairwise colony oceanographic niche distance')+ylab('Pairwise colony predictive performance')+
      annotate(geom = 'text', label = paste('Mantel r =', round(mn$statistic, 3), 'p =', round(mn$signif, 3)),
               x = Inf, y = Inf, hjust = 1, vjust = 1)+theme_bw()
    
    
  # niche plot
    enviro_std<-decostand(dat[,c(4:13)], method="standardize")
    enviro_rda<-rda(enviro_std, scale=T)
    #screeplot(enviro_rda)
    eig<-eigenvals(enviro_rda)
    enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1))
    enviro.sites.scores<-data.frame(enviro.sites.scores,dat[,c(1,2)])
    es1<-enviro.sites.scores%>%filter(forbin!= 'Core')%>%group_by(spcol)%>%slice(chull(PC1, PC2))%>%as.data.frame()
    es2<-enviro.sites.scores%>%filter(forbin!= 'Core')%>%
      group_by(spcol)%>%summarise(PC1=median(PC1), PC2=median(PC2))%>%as.data.frame()
    
    pniche<-ggplot()+
      geom_polygon(data=es1, aes(x=PC1, y=PC2, colour=spcol,fill=spcol), alpha=0.3)+
      geom_point(data=es2,aes(x=PC1, y=PC2, colour=spcol))+theme_bw()+
      scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100* eig[2]/sum(eig))))+
      scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)', 100* eig[1]/sum(eig))))
    
    multiplot<-p1+p2+pniche
      
      return(multiplot)
}
#### ~~~~ *** ~~~~ ####
#### ~~~~ Export mantel and niche plot ~~~~ ####

png('C:/seabirds/plots/brbo_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('BRBO', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/mabo_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('MABO', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/rfbo_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('RFBO', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/frbd_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('FRBD', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/trbd_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('TRBD', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/wtst_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('WTST', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/wtlg_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('WTLG', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/sote_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('SOTE', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/nodd_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('NODD', my.metric='AUC'))
dev.off()

png('C:/seabirds/plots/tern_mantel_niche.png',width = 12, height =4 , units ="in", res =300)
print(mklocada('TERN', my.metric='AUC'))
dev.off()


#### ~~~~ *** ~~~~ ####

#### ~~~~ Get env core/hull summary values ~~~~ ####
all_env<-NULL
for(i in c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
           'FRBD', 'TRBD', 'NODD', 'TERN'))
{
dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', i, '_kernhull_sample1.csv'))
dat$X<-NULL
dat$sp=i
all_env<-rbind(all_env, dat)
}

env_sum<-all_env%>%dplyr::select(-weight)%>%group_by(sp, spcol, forbin)%>%summarise_all(mean)
write.csv(env_sum, 'C:/seabirds/data/env_varib_mean_vals.csv', quote=F, row.names=F)

env_sum<-all_env%>%filter(forbin=='PsuedoA')%>%dplyr::select(-weight, -Longitude, -Latitude)%>%
  group_by(sp, spcol, forbin)%>%summarise(mn_sst=min(sst), mx_sst=max(sst),mn_chl=min(chl),
                                          mx_chl=max(chl),mn_mfr=min(mfr), mx_mfr=max(mfr),
                                          mn_pfr=min(pfr), mx_pfr=max(pfr),mn_chlsd=min(chl_sd),
                                          mx_chlsd=max(chl_sd),mn_sstsd=min(sst_sd), mx_sstsd=max(sst_sd),
                                          mn_mfrsd=min(mfr_sd), mx_mfrsd=max(mfr_sd),mn_pfrsd=min(pfr_sd),
                                          mx_pfrsd=max(pfr_sd),mn_bth=min(bth), mx_bth=max(bth),
                                          mn_slp=min(slp), mx_slp=max(slp))

plotdat<-all_env%>%filter(forbin=='PsuedoA')%>%dplyr::select(-weight, -Longitude, -Latitude, -forbin)%>%
  tidyr::gather('val', 'dat',-sp, -spcol)
# add GBR local data bounds: for each sp extract vals from GBR pred rasters that are within foraging range

for(i in c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG','FRBD', 'TRBD', 'NODD', 'TERN'))
{
  png(paste0('C:/seabirds/plots/env_varib_col_diffs_',i,'.png'),width = 11, height =6 , units ="in", res =300)
p1<-qplot(data=filter(plotdat, sp==i), y=dat, x=spcol, geom='boxplot')+
  facet_wrap(~val, scales='free', nrow=2)+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1),
                    axis.title.x=element_blank(),axis.title.y=element_blank())
print(p1)
dev.off()}
#### ~~~~ *** ~~~~ ####