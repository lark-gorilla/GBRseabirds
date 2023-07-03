# Code to accompany paper:

# Code takes outputs from random forest SDMs and transferability analyses and:
# calculates summary foraging radius values;
# conducts global Framework validation analyses;
# conducts foraging radius refinement on the GBR;
# calculates summary tables for paper
# generates non-GIS paper figures

# Mark Miller, mark.miller1@monash.edu, mark.gr.miller@gmail.com,

library(raster)
library(ggplot2)
library(ggspatial)
library(ggdendro)
library(dplyr)
library(tidyr)
library(ggplotify)
library(patchwork)
library(vegan)
library(sf)
library(viridis)
library(geodist)
library(smoothr)

setwd("C:/Users/mmil0049/Documents/sourced_data/MAML/C/seabirds")

####~~~~ read in data ~~~~####
aucz_out<-read.csv('data/mod_validation_vals_corrected.csv')
# auc tss cor
auc_temp<-filter(aucz_out, as.character(Resample)!=as.character(spcol))
auc_temp<-filter(auc_temp ,Resample!='EnsembleRaw')
auc_temp<-filter(auc_temp ,Resample!='EnsembleNrm')
auc_temp%>%group_by(sp)%>%summarise(cor(auc, TSS))
cor.test(auc_temp$auc, auc_temp$TSS)
qplot(data=auc_temp, x=auc, y=TSS)+facet_wrap(~sp)
# end cor test

matx_out<-read.csv('data/mod_clustering_vals.csv')
t_qual<-read.csv('data/tracking_trip_decisions.csv')
colz<-st_read('data/GIS/trackingID_colony_locs.shp')
spgroup_summ<-read.csv('data/sp_main_summary.csv')
gbr_cols<-read.csv('data/parks_gbr_colony_data.csv')
#### ~~~~ **** ~~~~ ####

####~~~~ species-colony dataset/foraging range suppl mat table ~~~~####

# formatting col
colz$sp<-do.call(c, lapply(strsplit(as.character(colz$ID), '_'), function(x)x[2]))
colz$coly<-do.call(c, lapply(strsplit(as.character(colz$ID), '_'), function(x)x[3]))
colz$spcol<-paste(colz$sp, colz$coly)
colz$Longitude<-st_coordinates(colz)[,1]
colz$Latitude<-st_coordinates(colz)[,2]
st_geometry(colz)<-NULL
colz<-colz%>%group_by(sp, coly)%>%summarise_all(first)
#png('C:/seabirds/outputs/maps/tracking_data_locs.png'),width = 8.3, height =8 , units ="in", res =300)
#ggplot()+geom_sf(data=land)+geom_sf(data=colz, shape=1, colour='red')+facet_wrap(~sp_group)+coord_sf(xlim = c(-180, 180), ylim = c(-36, 36), expand = F)+theme_bw()
#dev.off()

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

## breedstage fix
master<-read.csv('sourced_data/tracking_data/tracking_master.csv')
master$ID_tr<-paste(master$dataID, master$sp, master$colony, master$trackID, sep='_')
t_qual$ID_tr<-paste(t_qual$ID, t_qual$trackID, sep='_')

t_qual<-left_join(t_qual, master%>%group_by(ID_tr)%>%summarise_all(first)%>%dplyr::select(ID_tr,breedstage), by='ID_tr')
##
t_qual_unfilt<-t_qual

spcol_tab<-t_qual_unfilt%>%group_by(sp_group, sp, coly)%>%
  summarise(n_rawtracks=length(unique(trip_id)), col_long=first(Longitude), col_lat=first(Latitude))%>%as.data.frame()
# keep only manually filtered good trips
t_qual<-t_qual[t_qual$manual_keep=='Y',]
t_nused<-t_qual%>%group_by(sp_group, sp, coly)%>%
  summarise(n_usedtracks=length(unique(trip_id)))%>%as.data.frame()
spcol_tab<-data.frame(spcol_tab, n_usedtracks=t_nused$n_usedtracks)

# Returning trips as classified track2kba package
t_qual_ret<-t_qual[t_qual$complete=='complete trip',]
t_tripmetric<-t_qual_ret%>%filter(duration<296.1)%>%group_by(sp_group, sp, coly)%>%
  summarise(max_for=max(max_dist), mean_for=mean(max_dist), sd_for=sd(max_dist), breedstage.x=paste(unique(breedstage.x), collapse=' '),
            breedstage.y=paste(unique(breedstage.y), collapse=' '))%>%as.data.frame()
spcol_tab<-data.frame(spcol_tab, t_tripmetric[,c(4:8)])

# Returning trips are implicit as all GPS data recovered from colony, either from device on bird or near-colony download.
t_tripmetric2<-t_qual%>%filter(duration<296.1)%>%group_by(sp_group, sp, coly)%>%
  summarise(max_for_allret=max(max_dist), mean_for_allret=mean(max_dist), sd_for_allret=sd(max_dist))%>%as.data.frame()
spcol_tab<-data.frame(spcol_tab, t_tripmetric2[,c(4:6)])
# tidy
spcol_tab[spcol_tab$coly=='chick',]$coly<-'Rat' # edit to name
spcol_tab[spcol_tab$sp_group=='WTST' | spcol_tab$sp_group=='WTLG',]$sp<-'WTSH'
spcol_tab[spcol_tab$coly=='Adele' & spcol_tab$sp=='BRBO',]$max_for<-139.6  # manual edit
spcol_tab[spcol_tab$coly=='Adele' & spcol_tab$sp=='BRBO',]$max_for_allret<-139.6  # manual edit
spcol_tab[spcol_tab$coly=='Heron' & spcol_tab$sp_group=='WTLG',]$max_for_allret<-1150  # manual edit


#write.csv(spcol_tab, 'data/sp_col_summary_updatedJul23.csv', quote=F, row.names=F)

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

####~~~~ Paper main results table 3 ~~~~####

spcol_tab<-read.csv('data/sp_col_summary.csv')

for_rang<-spcol_tab%>%group_by(sp_group)%>%summarise(min.for=min(max_for_allret),
                                                mean.for=mean(max_for_allret),
                                                sd.for=sd(max_for_allret),
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


bind_out2<-bind_out%>%select(sp, mn_sumr, mx_sumr, mn_sumr_ts, mx_sumr_ts, auc_rank, tss_rank, auc_rank1, tss_rank1, min.for, mean.for, max.for)
#write.csv(bind_out2, 'C:/seabirds/data/sp_main_summary.csv', quote=F, row.names=F)

# extra multicol vs local model summary
my.aucz_temp<-aucz_out
my.aucz_temp<-left_join(my.aucz_temp, filter(my.aucz_temp, Resample=='MultiCol')%>%
                          select(sp, spcol, auc, TSS),by=c('sp', 'spcol'))
my.aucz_temp<-filter(my.aucz_temp, as.character(Resample)==as.character(spcol))

my.aucz_temp$self_mult_dclass_auc<-ifelse(abs(my.aucz_temp$auc.x-my.aucz_temp$auc.y)< 0.05, 'same',
                                          ifelse(my.aucz_temp$auc.x-my.aucz_temp$auc.y< 0, 'better', 'worse'))

my.aucz_temp$self_mult_dclass_tss<-ifelse(abs(my.aucz_temp$TSS.x-my.aucz_temp$TSS.y)< 0.05, 'same',
                                          ifelse(my.aucz_temp$TSS.x-my.aucz_temp$TSS.y< 0, 'better', 'worse'))

sumr1<-my.aucz_temp%>%group_by(sp, self_mult_dclass_auc)%>%summarise(nclass=n())

#### ~~~~ **** ~~~~ ####

#### ~~~~ Applying framework: refinement of foraging circles on the GBR ~~~~ ####
glob_auc<-data.frame(md_spgr=c('BRBO','MABO','RFBO','FRBD','TRBD','WTST','WTLG','SOTE','NODD','TERN'),
                     auc=c(0.55,0.53, 0.54, 0.61, 0.56, 0.58, 0.54, 0.48, 0.40, 0.82))

for_rad<-read_sf('data/GIS/foraging_radii.shp')

pred_list<-list.files('data/modelling/GBR_preds/selected_preds', full.names=T)
mod_pred<-stack(pred_list)
r_sp<-substr(pred_list[nchar(pred_list)==69], 53, 65) # selects Multi-colony models
r_sp<-r_sp[c(9,2,8,6,5,10,3,4,1,7)] # order by range for memory troubleshooting run

collect_polys<-NULL
for( i in r_sp)
{
  spkey=substr(i, 1, 4)
  col<-filter(for_rad, md_spgr==spkey & rd_clss %in% c('mea', 'obs'))
  
  for(j in unique(col$site_nm))
  {
    sp_ras1<-subset(mod_pred, i)
    col_sp<-filter(col, site_nm==j)
    
    if('obs' %in% col_sp$rd_clss)
    {
      col_sp<-filter(col_sp, rd_clss=='obs')
      if(spkey=='BRBO' & col_sp$dsgntn_n=='Raine Island, Moulter and MacLennan cays KBA'){
        sp_ras1<-subset(mod_pred, 'BRBO_Raine')}
      if(spkey=='BRBO' & col_sp$dsgntn_n=='Swain Reefs KBA'){
        sp_ras1<-subset(mod_pred, 'BRBO_Swains')}
      if(spkey=='MABO' & col_sp$dsgntn_n=='Swain Reefs KBA'){
        sp_ras1<-subset(mod_pred, 'MABO_Swains')}
      if(spkey=='WTST' & col_sp$dsgntn_n=='Capricornia Cays KBA'){
        sp_ras1<-subset(mod_pred, 'WTST_Heron')}
      if(spkey=='WTLG' & col_sp$dsgntn_n=='Capricornia Cays KBA'){
        sp_ras1<-subset(mod_pred, 'WTLG_Heron')}
      if(spkey=='NODD' & col_sp$dsgntn_n=='Capricornia Cays KBA'){
        sp_ras1<-subset(mod_pred, 'NODD_Heron')}
      # Not for WTLG Heron as good already
    }
    
    sp_ras<-crop(sp_ras1, extent(col_sp)) # drop size
    
    sum_rad<-mask(sp_ras, as(col_sp, 'Spatial'), updatevalue=999, updateNA=T)
    
    sr2<-sum_rad
    sr2[sr2==999]<-NA
    rm(sum_rad)
    rm(sp_ras)
    rm(sp_ras1)
    
    # lookup auc global
    auc_val<-glob_auc[glob_auc$md_spgr==col_sp$md_spgr,]$auc
    
    if(col_sp$md_spgr=='BRBO' & col_sp$dsgntn_n=='Raine Island, Moulter and MacLennan cays KBA'){
      auc_val<-0.66}
    if(col_sp$md_spgr=='BRBO' & col_sp$dsgntn_n=='Swain Reefs KBA'){auc_val<-0.65}
    if(col_sp$md_spgr=='MABO' & col_sp$dsgntn_n=='Swain Reefs KBA'){auc_val<-0.64}
    if(col_sp$md_spgr=='WTST' & col_sp$dsgntn_n=='Capricornia Cays KBA'){auc_val<-0.74}
    if(col_sp$md_spgr=='WTLG' & col_sp$dsgntn_n=='Capricornia Cays KBA'){auc_val<-0.64}
    if(col_sp$md_spgr=='NODD' & col_sp$dsgntn_n=='Capricornia Cays KBA'){auc_val<-0.66}
    
    #seq1<-seq(auc_val, 0.5+auc_val, 0.1)
    #seq1<-seq1[seq1< 1 & seq1>=0.51] 
    #seq1<-seq1[seq1!=auc_val]
    
    auc_trials<-c(auc_val, seq(0.5, 0.9, 0.02))
    for(k in auc_trials)
    {
      #scale auc to quantile cutoff: 0.5=0 (ie radius), 0.9=0.9 ('Excellent' prediction gives 10% core)
      auc_cut<-((k-0.5)/(0.9-0.5))*(0.9-0)+0
      if(auc_cut>0.9){auc_cut<-0.9} # anything above 0.9AUC gets 10% core also
      if(auc_cut<0){auc_cut=0} # catch <0.5AUC species and set to radius
      
      q2 <- quantile(sr2, auc_cut)
      ## ** ##
      hots<-reclassify(sr2, c(-Inf, q2, NA, q2, Inf, 1), right=F)
      s1<-st_as_sf(rasterToPolygons(hots, dissolve = T)) 
      rm(hots)
      names(s1)[1]<-'mod'
      s2<-drop_crumbs(s1, threshold=16000000)
      rm(s1)
      s3<-fill_holes(s2, threshold=54000000)
      rm(s2)
      s3 <- smoothr::smooth(s3, method = "ksmooth", smoothness = 4, n=1) # n=1 stops over-densifying
      s3<-drop_crumbs(s3, threshold=16000000)

      # rbind sf object
      s3$dsgntn_n<-col_sp$dsgntn_n[1]
      s3$site_nm<-col_sp$site_nm[1]
      s3$md_spgr<-col_sp$md_spgr[1]
      s3$species<-col_sp$species[1]
      s3$auc_type<-'sim'
      if(which(k==auc_trials)==1){s3$auc_type<-'obs'}
      s3$auc<-k
      s3$area_km2<-as.numeric(st_area(s3)/1000000)
      s3$mod<-NULL
      #cp2<-st_read('C:/seabirds/temp/auc_cores_temp.shp')
      #cp2<-rbind(cp2, s3)
      #st_write(cp2, 'C:/seabirds/temp/auc_cores_temp.shp', delete_dsn = T)
      rm(cp2)
      collect_polys <- rbind(collect_polys, s3)
      rm(s3)
      
    }
  }
print(i)
}

# write polys
#st_write(collect_polys, 'C:/seabirds/data/GIS/col_radii_auc_core_smooth_sim0.2auc.shp', delete_dsn=T)
#### ~~~~ **** ~~~~ ####

#### ~~~~ Framework validation: known foraging area inclusion in refined foraging circles ~~~~ ####
# rad prediction data
rad_pred<-read.csv('data/radius_LGOCV_predictions_2km_dynamic.csv')
names(rad_pred)[1]<-'ID'


templ<-raster('sourced_data/oceano_modelready/extraction_template_2km_chl.tif')


# kernels
kernz<-rbind(st_read('data/GIS/NODDkernhull.shp'),
             st_read('data/GIS/SOTEkernhull.shp'),
             st_read('data/GIS/TERNkernhull.shp'),
             st_read('data/GIS/BRBOkernhull.shp'),
             st_read('data/GIS/MABOkernhull.shp'),
             st_read('data/GIS/RFBOkernhull.shp'),
             st_read('data/GIS/FRBDkernhull.shp'),
             st_read('data/GIS/TRBDkernhull.shp'))
kernz$spcol<-as.character(kernz$spcol)
kernz[kernz$spcol=='SOTE chick',]$spcol<-'SOTE Rat'

w1<-st_read('data/GIS/WTSTkernhull.shp')
w1$spcol<-gsub( 'WTSH', 'WTST', w1$spcol)
kernz<-rbind(kernz, w1)
rm(w1)
w1<-st_read('data/GIS/WTLGkernhull.shp')
w1$spcol<-gsub( 'WTSH', 'WTLG', w1$spcol)
kernz<-rbind(kernz, w1)

#### auc vals
aucz_out<-read.csv('data/mod_validation_vals_corrected.csv')
aucz_out$spcol<-as.character(aucz_out$spcol)
aucz_out[aucz_out$spcol=='Lord Howe',]$spcol<-'LHI'
aucz_out[aucz_out$sp=='TRBD' & aucz_out$spcol=='Pena Blanca',]$spcol<-'Peña blanca'

# read for rads
for_rad<-read_sf('data/GIS/global_mean_foraging_radii.shp')
for_rad$spcol<-paste(for_rad$sp, for_rad$coly) # wedgies

cover_tab<-NULL
for(i in unique(rad_pred$ID))
{
  if(i=="RFBO Christmas"){next}
  if(i=="RFBO Isabel"){next}
  if(i=="WTLG Aride"){next}
  i=as.character(i)
  dat1<-rad_pred[rad_pred$ID==i,]
  
  #lookup up kern
  for_kern<-filter(kernz, spcol==i & PA==1)
  for_kern<-as(for_kern, 'Spatial')
 
  spdf<-SpatialPointsDataFrame(SpatialPoints(dat1[c('x', 'y')],
                proj4string = CRS(projection(templ))), data=dat1['pred'])
  
  # make combined extent
  comb_bbox<-st_union(st_as_sfc(st_bbox(spdf)),
              st_as_sfc(st_bbox(for_kern)))
  
  sp_templ<-crop(templ, extent(as(comb_bbox, 'Spatial')))
  
  kern_ras<-rasterize(for_kern, sp_templ, field=1,background=0)
  
  
  # radius raster
  rad_ras<-rasterize(spdf, sp_templ, field=1,background=0)
  
  # radius dist raster - option to weight over distance to col cost surface - NOT USED
  #dist1<-rasterize(as(st_centroid(for_rad[for_rad$spcol==i,]), 'Spatial'),
  #                 sp_templ, field=1,background=NA)
  
  #dist1<-distance(dist1)
  #dist1<-mask(dist1, rad_ras, maskvalue=0) # mask out land
  #dist1<-calc(dist1, fun=function(x){abs(x-max(values(dist1), na.rm=T))})
  
  # no calc refined rad raster
  r1<-rasterize(spdf, sp_templ, field='pred') # set background to NA
  
  plot(r1, main=i)
  plot(for_kern, add=T)
  # filter auc by sp first then to spcol
  sp_filt<-substr(paste(i), 1, 4)
  if(sp_filt %in% c('RTTB', 'RBTB')){sp_filt<-'TRBD'}
  if(sp_filt %in% c('CRTE', 'CATE', 'ROTE')){sp_filt<-'TERN'}
  if(sp_filt %in% c('MAFR', 'LEFR', 'GRFR')){sp_filt<-'FRBD'}
  if(sp_filt %in% c('BLNO', 'BRNO', 'LENO')){sp_filt<-'NODD'}
  
  sp_auc<-aucz_out[aucz_out$sp==sp_filt,]
  i_lkup<-substr(i,6, nchar(i))
  if(sp_filt=='TERN'|sp_filt=='FRBD'){i_lkup<-i}
  
  auc_lookup<-sp_auc[sp_auc$Resample=='MultiCol' & sp_auc$spcol==i_lkup,]$auc
  
  seq1<-c(auc_lookup, 0.5,0.51,0.56,0.61,0.66,0.71,0.76,
          0.81,0.86, 0.91)# add 0.05 intervals for full run
  # remember first val in output table is mod supported
  for(k in seq1)
  {

    auc_cut<-((k-0.5)/(0.9-0.5))*(0.9-0)+0
    if(auc_cut>0.9){auc_cut<-0.9} # anything above 0.9AUC gets 10% core also
    if(auc_cut<0){auc_cut=0} # catch <0.5AUC species and set to radius
    
    if (k<0.51){ref_ras=rad_ras; dist_ras=rad_ras}else
    { 
    #refining model pred
    q2 <- quantile(r1, auc_cut)
    hots<-reclassify(r1, c(-Inf, q2, NA, q2, Inf, 1), right=F)
    
    # same refinement process in smoothing and dropping crumbs
    s1<-st_as_sf(rasterToPolygons(hots, dissolve = T)) 
    rm(hots)
    names(s1)[1]<-'mod'
    s2<-drop_crumbs(s1, threshold=16000000)
    rm(s1)
    s3<-fill_holes(s2, threshold=54000000)
    rm(s2)
    s3 <- smoothr::smooth(s3, method = "ksmooth", smoothness = 4, n=1) # n=1 stops over-densifying
    s3<-drop_crumbs(s3, threshold=16000000)
    # and rasterize
    ref_ras<-rasterize(as(s3, 'Spatial'), sp_templ, field=1,background=0)
    
    ## refining dist mod
    #t1<-as.data.frame(table(values(dist1)))
    #t1<-t1[order(t1$Var1, decreasing = T),]
    #t1$cum_freq<-cumsum(t1$Freq)
    #t1$diff=abs(t1$cum_freq-table(values(ref_ras))['1'])
    
    #qval<-as.numeric(as.character(t1[which.min(t1$diff),]$Var1)) 
 
    #hots<-reclassify(dist1, c(-Inf, qval, NA, qval, Inf, 1), right=F)
    #s3<-st_as_sf(rasterToPolygons(hots, dissolve = T)) # no need to smooth or drop holes in radius
    #rm(hots)
    #dist_ras<-rasterize(as(s3, 'Spatial'), sp_templ, field=1,background=0)
    }
  
  #mask out land
  kern_ras<-mask(kern_ras, sp_templ)
  rad_ras<-mask(rad_ras, sp_templ)
  ref_ras<-mask(ref_ras, sp_templ)
  #dist_ras<-mask(dist_ras, r1)
  
  #calc percentage overlap
  npix_kern<-table(values(kern_ras))['1']
  npix_kern_inrad<-table(values(kern_ras+rad_ras))['2']
  npix_kern_inref<-table(values(kern_ras+ref_ras))['2']
  #npix_kern_indist<-table(values(kern_ras+dist_ras))['2']
  
  plot(stack(kern_ras+rad_ras, kern_ras+ref_ras), main=i)
  
  d1<-data.frame(sp=sp_filt, spcol=i, auc=k, perc_cut=auc_cut,
                 percfor_inrad=round(npix_kern_inrad/npix_kern*100),
                 percfor_inref=round(npix_kern_inref/npix_kern*100),
                 npixkern=npix_kern, npixrad=table(values(rad_ras))['1'],
                 npixref=table(values(ref_ras))['1'])
  
  cover_tab<-rbind(cover_tab, d1)
  } #close k loop 
  print(i)
}

#write.csv(cover_tab, 'data/refined_radii_core_inclusion_2km_dynamic.csv', quote=F, row.names=F)
#### ~~~~ **** ~~~~ ####

#### ~~~~ Local GBR tracking known foraging area inclusion in refined foraging circles ~~~~ ####

for_rad<-read_sf('/data/GIS/foraging_radii.shp')
for_rad<-filter(for_rad, rd_clss=='obs')
for_rad<-rbind(
  filter(for_rad, md_spgr=='BRBO' & site_nm=='Raine Island'),
  filter(for_rad, md_spgr=='BRBO' & site_nm=='Price Cay'),
  filter(for_rad, md_spgr=='MABO' & site_nm=='Price Cay'),
  filter(for_rad, md_spgr=='NODD' & site_nm=='Heron Island'),
  filter(for_rad, md_spgr=='WTST' & site_nm=='Heron Island'),
  filter(for_rad, md_spgr=='WTLG' & site_nm=='Heron Island'))
  

pred_list<-list.files('data/modelling/GBR_preds/selected_preds', full.names=T)
mod_pred<-stack(pred_list)

cover_tab<-NULL
for( i in 1:nrow(for_rad))
{
  col<-for_rad[i,]
  spkey<-col$md_spgr
  
      if(spkey=='BRBO' & col$site_nm=='Raine Island'){
        sp_ras1<-subset(mod_pred, 'BRBO_Raine');auc_val<-0.66;spc<-'BRBO Raine'}
      if(spkey=='BRBO' & col$site_nm=='Price Cay'){
        sp_ras1<-subset(mod_pred, 'BRBO_Swains'); auc_val<-0.65;spc<-'BRBO Swains'}
      if(spkey=='MABO'){
        sp_ras1<-subset(mod_pred, 'MABO_Swains');auc_val<-0.64;spc<-'MABO Swains'}
      if(spkey=='WTST'){
        sp_ras1<-subset(mod_pred, 'WTST_Heron');auc_val<-0.74;spc<-'WTST Heron'}
      if(spkey=='WTLG'){
        sp_ras1<-subset(mod_pred, 'WTLG_Heron');auc_val<-c(0.64, 0.9);spc<-'WTLG Heron'}
      if(spkey=='NODD'){
        sp_ras1<-subset(mod_pred, 'NODD_Heron');auc_val<-0.66;spc<-'BLNO Heron'}
    
    #clip pred inside radius  
    sp_ras<-crop(sp_ras1, extent(col)) # drop size
    sum_rad<-mask(sp_ras, as(col, 'Spatial'), updatevalue=999, updateNA=T)
    sr2<-sum_rad
    sr2[sr2==999]<-NA
    rm(sum_rad)
    rm(sp_ras1)
    
    #lookup up kern
    for_kern<-filter(kernz, spcol==spc & PA==1)
    for_kern<-as(for_kern, 'Spatial')
    kern_ras<-rasterize(for_kern, sp_ras, field=1,background=0)
    
    for(k in auc_val)
      {
        
        auc_cut<-((k-0.5)/(0.9-0.5))*(0.9-0)+0
        if(auc_cut>0.9){auc_cut<-0.9} # anything above 0.9AUC gets 10% core also
        if(auc_cut<0){auc_cut=0} # catch <0.5AUC species and set to radius
        
          #refining model pred
          q2 <- quantile(sr2, auc_cut)
          hots<-reclassify(sr2, c(-Inf, q2, NA, q2, Inf, 1), right=F)
          
          # same refinement process in smoothing and dropping crumbs
          s1<-st_as_sf(rasterToPolygons(hots, dissolve = T)) 
          rm(hots)
          names(s1)[1]<-'mod'
          s2<-drop_crumbs(s1, threshold=16000000)
          rm(s1)
          s3<-fill_holes(s2, threshold=54000000)
          rm(s2)
          s3 <- smoothr::smooth(s3, method = "ksmooth", smoothness = 4, n=1) # n=1 stops over-densifying
          s3<-drop_crumbs(s3, threshold=16000000)
          # and rasterize
          ref_ras<-rasterize(as(s3, 'Spatial'), sp_ras, field=1,background=0)

        #mask out land
        kern_ras<-mask(kern_ras, sp_ras)
        ref_ras<-mask(ref_ras, sp_ras)
        #dist_ras<-mask(dist_ras, r1)
        
        #calc percentage overlap
        npix_kern<-table(values(kern_ras))['1']
        npix_kern_inref<-table(values(kern_ras+ref_ras))['2']

        plot(kern_ras+ref_ras, main=i)
        
        d1<-data.frame(spcol=spc, auc=k, perc_cut=auc_cut,
                       percfor_inref=round(npix_kern_inref/npix_kern*100),
                       npixkern=npix_kern,
                       npixref=table(values(ref_ras))['1'])
        
        cover_tab<-rbind(cover_tab, d1)
    }
  print(i)
}
#write.csv(cover_tab, 'data/gbr_local_core_foraging_area_inclusion.csv', quote=F, row.names=F)

#### ~~~~ **** ~~~~ ####

#### ~~~~ Make framework validation plots and tables (Fig 4, Table 4, suppl fig~~~~ ####

#summary and plots
inclu<-read.csv('data/refined_radii_core_inclusion_2km_dynamic.csv')
inclu[is.na(inclu$percfor_inref),]$percfor_inref<-0 # replace NA with 0
inclu$ref_loss=inclu$percfor_inrad-inclu$percfor_inref # % loss compared to radius (ignores core foraging outside radius)
# group into AUC accuracy classes
inclu<-inclu%>%group_by(spcol)%>%mutate(min_auc=first(auc), auc_bin=cut(min_auc, breaks=c(-Inf, 0.5, 0.6, 0.7, 0.8, 0.9, Inf),
                        labels=c('radius', 'v poor', 'poor', 'moderate', 'good', 'excellent')))

#select observed data
inclu_obs<-inclu%>%group_by(spcol)%>%summarise_all(first)%>%as.data.frame()
#select simulated data
inclu_sim<-inclu%>%group_by(spcol)%>%slice(2:n())%>%as.data.frame()# get bottom 6 rows

ggplot(data=inclu_sim,
      aes(x=auc, y=percfor_inref/100, colour=auc_bin))+
  geom_smooth(method = "glm", method.args = list(family = "binomial"))


m_rad<-glmer(cbind(percfor_inref, 100-percfor_inref)~perc_cut+perc_cut:min_auc+(perc_cut|spcol), 
             data=inclu_sim, family='binomial', 
             control = glmerControl(optimizer ="bobyqa"))
print(sum((resid( m_rad, type="pearson")^2))/df.residual( m_rad))
plot(m_rad)
Anova(m_rad)

# Interpret log odds
#https://stats.stackexchange.com/questions/377515/computation-and-interpretation-of-odds-ratio-with-continuous-variables-with-inte
#https://www.reddit.com/r/statistics/comments/4begf8/question_about_interpreting_negative_coefficients/d18gjog/
summary(m_rad)$coefficients
exp(abs(-0.14768892+(0.12598312*c(0.5, 0.6,0.7,0.8,0.9))))


new_dat<-expand.grid(percfor_inref=1, perc_cut=seq(0, 0.9, 0.02), 
                     min_auc=c(0.5,0.6,0.7,0.8,0.9))
new_dat$pred<-predict(m_rad, new_dat, type='link', re.form=NA)

predmat <- model.matrix(terms(m_rad), data=new_dat) # need to put terms arguement not just model, RE carried over otherwise?
vcv <- vcov(m_rad)
## then calculate the standard errors
semod <-  sqrt(diag(predmat%*%vcv%*%t(predmat))) #M: creates matrix and takes the diagonal
# then we can get the confidence intervals @ 95% confidence level
new_dat$ucl <- new_dat$pred + semod*1.96
new_dat$lcl <- new_dat$pred - semod*1.96

# get observed data mn+sd 
obs_mn<-inclu_obs%>%filter(auc_bin!='radius')%>%group_by(auc_bin)%>%
  summarise(mn_r=mean(percfor_inref),sd_r=sd(percfor_inref))
inc_notrad<-inclu_obs[inclu_obs$auc_bin!='radius',]

#obs_mod<-glm(cbind(percfor_inref, 100-percfor_inref)~perc_cut, 
 #            data=inc_notrad[-c(34, 45),], family='binomial')
obs_mod<-gam(cbind(percfor_inref, 100-percfor_inref)~s(perc_cut, k=5), 
             data=inc_notrad[-c(34, 45),], family='binomial')

print(sum((resid( obs_mod, type="pearson")^2))/df.residual( obs_mod))

new_dat2<-expand.grid(percfor_inref=1, perc_cut=seq(0, 0.9, 0.02))
new_dat2<-cbind(new_dat2, predict(obs_mod, new_dat2, type='link', se.fit=T))

# for main
p1<-ggplot(data=new_dat, aes(x=perc_cut))+
     geom_point(data=inc_notrad, aes(y=percfor_inref/100), shape=1, alpha=0.5)+
     geom_ribbon(aes(ymin=plogis(lcl), ymax=plogis(ucl), fill=factor(min_auc)),alpha=0.3, colour=NA)+
     geom_line(aes(y=plogis(pred),colour=factor(min_auc)), size=1)+
    geom_ribbon(data=new_dat2, aes(ymin=plogis(fit- se.fit*1.96), ymax=plogis(fit+ se.fit*1.96)),
      alpha=0.3, fill=NA, colour='black', linetype='dashed')+
      geom_line(data=new_dat2, aes(y=plogis(fit)))+
     scale_x_continuous(breaks=c(0, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))+
     ylab('Probability of known foraging areas included in refined foraging circle')+
     xlab('Predicted habitat suitability threshold percentile')+
     guides(colour=guide_legend(title='Model transferability (AUC)'),
    fill=guide_legend(title='Model transferability (AUC)'))+theme_bw()+
     theme(legend.position = c(0.22, 0.17), legend.background = element_rect(fill=NA, 
     size=0.5, linetype="solid"))+
  annotate("text", x = 0.45, y = 0.025, label = "y=4.6006+(-14.7689*x)+(12.5983*x*AUC)")

#ggsave(plot=p1, width =5 , height =5, units='in',
#       filename='plots/all_sp_rad_refine_gam.png')

# make for each species seperately 
plot_dat<-NULL
coefz<-NULL
for(i in unique(inclu_sim$sp))
{
   print(i)
  sp_inclu<-inclu_sim[inclu_sim$sp==i,]
  sp_inclu$spcol<-factor(sp_inclu$spcol)
  sp_inclu$auc_bin<-factor(sp_inclu$auc_bin)
  
  # allow RE level intercept for these models
  if(i=='NODD'){
  m_sp<-glmer(cbind(percfor_inref, 100-percfor_inref)~perc_cut+(perc_cut|spcol), 
                data=sp_inclu[-12,], family='binomial', 
                control = glmerControl(optimizer ="bobyqa"))# NODD only has radius auc_bin 
    
  }else{
  m_sp<-glmer(cbind(percfor_inref, 100-percfor_inref)~perc_cut+perc_cut:min_auc+(perc_cut|spcol), 
               data=sp_inclu, family='binomial', 
               control = glmerControl(optimizer ="bobyqa"))}
  
  print(sum((resid( m_sp, type="pearson")^2))/df.residual( m_sp))
  print(plot(m_sp))
  readline('')
  new_dat<-expand.grid(sp=i,perc_cut=seq(0, 0.9, 0.02), 
                       min_auc=seq(min(sp_inclu$min_auc), max(sp_inclu$min_auc), 0.02))
  new_dat$pred<-predict(m_sp, new_dat, type='response', re.form=NA)
  
  # fit glm instead for sp wih singularity due to few REs
  if( i=='TERN'){
    m_sp<-glm(cbind(percfor_inref, 100-percfor_inref)~perc_cut+perc_cut:min_auc, 
                  data=sp_inclu, family='binomial') 
    print(sum((resid( m_sp, type="pearson")^2))/df.residual( m_sp))
    new_dat$pred<-predict(m_sp, new_dat, type='response')}
  if(i=='SOTE'){
    m_sp<-glm(cbind(percfor_inref, 100-percfor_inref)~perc_cut+perc_cut:min_auc, 
              data=sp_inclu, family='quasibinomial') 
    new_dat$pred<-predict(m_sp, new_dat, type='response')}
  
  print(ggplot(data=new_dat, aes(x=perc_cut, y=pred, group=min_auc,colour=min_auc))+
          geom_line()+
          geom_point(data=sp_inclu, aes(x=perc_cut, y=percfor_inref/100), alpha=0.5))
  
  readline('')
  plot_dat<-rbind(plot_dat, new_dat)
  
  c1<- data.frame(sp=i,coefs=summary(m_sp)$coefficients[,1],
                  sd=summary(m_sp)$coefficients[,2], 
                  var=row.names(summary(m_sp)$coefficients))
  coefz<-rbind(coefz, c1)
}

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
binreg_supl<-ggplot(data=plot_dat[plot_dat$min_auc>0.4,], aes(x=perc_cut, y=pred, colour=min_auc))+
        geom_line(aes(group=min_auc))+facet_wrap(~sp)+
  scale_colour_gradientn(colours=myPalette(10),breaks=c(0.5,0.6,0.7,0.8,0.9),
                         limits=c(0.4,1))+
  ylab('Probability of core foraging areas included in refined radius')+
  xlab('Predicted habitat suitability threshold percentile')+
  guides(colour=guide_legend(title='Model transferability (AUC)'))

#ggsave(binreg_supl,  width =8 , height =8, units='in',
#       filename='C:/seabirds/plots/bin_regress_sp_suppl.png')

# Interpret log odds
#https://stats.stackexchange.com/questions/377515/computation-and-interpretation-of-odds-ratio-with-continuous-variables-with-inte
#https://www.reddit.com/r/statistics/comments/4begf8/question_about_interpreting_negative_coefficients/d18gjog/
cof1<-data.frame(sp=coefz[coefz$var=='(Intercept)',]$sp, 
           perc_cut=coefz[coefz$var=='perc_cut',]$coefs, 
           perc_cut_min_auc=c(0, coefz[coefz$var=='perc_cut:min_auc',]$coefs))

cof1$oddsr_0.5AUC<-exp(abs((cof1$perc_cut/100)+((cof1$perc_cut_min_auc/100)*0.5)))
cof1$oddsr_0.7AUC<-exp(abs((cof1$perc_cut/100)+((cof1$perc_cut_min_auc/100)*0.7)))
cof1$oddsr_0.8AUC<-exp(abs((cof1$perc_cut/100)+((cof1$perc_cut_min_auc/100)*0.8))) 
# format coefz
coefz$text=paste0(round(coefz$coefs, 2), 
       '±',round(coefz$sd, 2))
m1<-tidyr::spread(coefz[,c(1,4,5)], key=var,value=text,fill = '')

#write.csv(m1, 'data/bin_regress_sp_coefs_cont_int.csv', quote=F, row.names=F)

# range of percentiles used for training
inclu_obs%>%group_by(sp)%>%summarise(minauc=min(min_auc), 
maxauc=max(min_auc), minp=min(perc_cut), maxp=max(perc_cut))

# get inclusion of known areas by radii

intz<-coefz[coefz$var=='(Intercept)',]
rad_inc<-data.frame(sp=intz$sp, pred=plogis(intz$coefs), 
      uci=plogis(intz$coefs+intz$sd*1.96), lci=plogis(intz$coefs-intz$sd*1.96))
rad_inc$text<-paste0(round(rad_inc$pred, 2), 
                     ' (',round(rad_inc$lci, 2),',',round(rad_inc$uci, 2), ')')

qplot(data=inclu_obs, x=auc, y=ref_loss)+facet_wrap(~sp)

p1<-ggplot(data=inclu_obs, aes(x=auc_bin, y=ref_loss))+
  geom_boxplot()+geom_jitter(height=0, width = 0.1, shape=1, colour='red', alpha=0.6)+
facet_wrap(~sp, scales='free')+ylab('Percentage of core foraging area excluded')+
  xlab('Radius refinement grouped by Multi-Colony model transferability (AUC)')

#ggsave(p1,  width =8 , height =8, units='in',
#       filename='plots/for_refine_rad_inc.png')

#### ~~~~ **** ~~~~ ####

#### ~~~~ GBR refined foraging circles for paper: Fig 5 and Table 5, and suppl table ~~~~ ####

# read in SIMPLIFIED auc-calced core areas
corez<-st_read('data/GIS/col_radii_auc_core_smooth_sim0.2auc_simp.shp')
corez<-corez%>%filter(!duplicated(.)) #remove 'obs' rows that are also sim 'rows'

#! REMEBER to RUN once !#
#! remove Wedgie Mudjimba colony
#corez<-corez[corez$dsgntn_!='Mudjimba Island',]
#! add perc_cut used
#corez$perc_cut<-((corez$auc-0.5)/(0.9-0.5))*(0.9-0)+0
#corez[corez$perc_cut<0,]$perc_cut<-0 
#! And simplify polys, recalc area and write out
#c_simpl<-corez%>%st_simplify(preserveTopology=T,dTolerance=0.007)
#c_simpl$are_km2<-as.numeric(st_area(c_simpl)/1000000)
#c_simpl$area_km2<-NULL
#st_write(c_simpl,'C:/seabirds/data/GIS/col_radii_auc_core_smooth_sim0.2auc_simp.shp', delete_dsn=T)

# summarise cores
corez_nosf<-corez%>%as.data.frame()
corez_nosf$geometry<-NULL

# add mod type
corez_nosf$mod<-'global'
corez_nosf[corez_nosf$md_spgr=='BRBO' & corez_nosf$dsgntn_n=='Raine Island, Moulter and MacLennan cays KBA',]$mod<-'local'
corez_nosf[corez_nosf$md_spgr=='BRBO' & corez_nosf$dsgntn_n=='Swain Reefs KBA',]$mod<-'local2'
corez_nosf[corez_nosf$md_spgr=='MABO' & corez_nosf$dsgntn_n=='Swain Reefs KBA',]$mod<-'local'
corez_nosf[corez_nosf$md_spgr=='WTST' & corez_nosf$dsgntn_n=='Capricornia Cays KBA',]$mod<-'local'
corez_nosf[corez_nosf$md_spgr=='NODD' & corez_nosf$dsgntn_n=='Capricornia Cays KBA',]$mod<-'local'
corez_nosf[corez_nosf$md_spgr=='WTLG' & corez_nosf$dsgntn_n=='Capricornia Cays KBA',]$mod<-'local'

# add site-level point colours
corez_nosf$st_c_ty<-as.character(corez_nosf$auc_type)
# first auc val below 100000 for some sp
lookup1<-corez_nosf%>%filter(md_spgr%in%unique(corez[corez$are_km2>100000,]$md_spgr))%>%
filter(are_km2<100000)%>%group_by(md_spgr, site_nm)%>%summarise_all(first)
# manually add wtlg cap cays that are above 100k km2 eve at 10% core
lookup1<-bind_rows(lookup1, filter(corez_nosf, md_spgr=='WTLG' & dsgntn_n=='Capricornia Cays KBA' & auc==0.9))
corez_nosf[paste(corez_nosf$md_spgr, corez_nosf$site_nm, corez_nosf$auc) %in%
           paste(lookup1$md_spgr, lookup1$site_nm, lookup1$auc),]$st_c_ty<-'first'
corez_nosf$st_c_ty<-ifelse(corez_nosf$auc_type=='obs', 'obs', corez_nosf$st_c_ty)

#! make main text summary table !#
# using md_spgr dissolved polys

## dissolve
corez$st_c_ty<-corez_nosf$st_c_ty #attrib 
rad_diss<-corez%>%filter(auc==0.5)%>%group_by(md_spgr)%>%summarise(geometry=st_union(geometry))
obs_diss<-corez%>%filter(auc_type=='obs')%>%group_by(md_spgr)%>%summarise(geometry=st_union(geometry))
conf_diss<-corez%>%filter(st_c_ty!='sim')%>%group_by(md_spgr, site_nm)%>%
  filter(st_c_ty==last(st_c_ty))%>%ungroup()%>%group_by(md_spgr)%>%summarise(geometry=st_union(geometry))

rad_diss$are_km2<-as.numeric(st_area(rad_diss)/1000000)
obs_diss$are_km2<-as.numeric(st_area(obs_diss)/1000000)
conf_diss$are_km2<-as.numeric(st_area(conf_diss)/1000000)

tab1<-data.frame(filter(corez_nosf, auc==0.5 )%>%group_by(md_spgr)%>%
  summarise(n_site=length(unique(dsgntn_n)), n_col= length(unique(site_nm))), 
            rad_area=rad_diss$are_km2/1000,obs_area=obs_diss$are_km2/1000,
            conf_area=conf_diss$are_km2/1000)

# print all sp-dissolved colum totals
(as.numeric(rad_diss%>%st_union%>%st_area)/1000000)/1000 #2941
(as.numeric(obs_diss%>%st_buffer(0)%>%st_union%>%st_area)/1000000)/1000 #2744
(as.numeric(conf_diss%>%st_buffer(0)%>%st_union%>%st_area)/1000000)/1000 #1115

# Attribute confidence to GBR predictions: 1) make suppl table per colony; 2) summarise for tab1

# read in coefs to look up
coefz<-read.csv('C:/seabirds/data/bin_regress_sp_coefs_cont_int.csv')
coefz$int<-as.numeric(unlist(lapply(strsplit(as.character(coefz$`X.Intercept.`), '±'), function(x){x[1]})))
coefz$pc<-as.numeric(unlist(lapply(strsplit(as.character(coefz$perc_cut), '±'), function(x){x[1]})))
coefz$pc_auc<-as.numeric(unlist(lapply(strsplit(as.character(coefz$perc_cut.min_auc), '±'), function(x){x[1]})))
coefz$pc_auc[1]<-0

# confidence of observed vs observed + area_target (some)
nosf_obs<-corez_nosf%>%filter(auc_type=='obs')
nosf_conf<-corez_nosf%>%filter(st_c_ty!='sim')%>%group_by(md_spgr, site_nm)%>%
  filter(st_c_ty==last(st_c_ty))%>%ungroup()

# obs % core foraging areas predicted in radii refined using auc-perc specified refinement
nosf_obs<-left_join(nosf_obs, coefz[,c(1,5,6,7)], by=c('md_spgr'='sp'))
nosf_obs$prob_cfor_inc<-plogis(nosf_obs$int+ (nosf_obs$pc*nosf_obs$perc_cut)+
                             (nosf_obs$pc_auc*nosf_obs$perc_cut*nosf_obs$auc)) 

# obs + confidence % core foraging areas predicted in radii refined using auc-perc specified refinement
# first we look up percentiles used from confidence 
nosf_conf<-left_join(nosf_conf, coefz[,c(1,5,6,7)], by=c('md_spgr'='sp')) 
# then we get observed model auc to show the original model transferability
nosf_conf$auc_obs<-nosf_obs$auc
nosf_conf$prob_cfor_inc<-plogis(nosf_conf$int+ (nosf_conf$pc*nosf_conf$perc_cut)+
                                 (nosf_conf$pc_auc*nosf_conf$perc_cut*nosf_conf$auc_obs)) 

# make supplement table of all colonies
supl_inc<-nosf_obs
supl_inc$perc_cut_areaT<-nosf_conf$perc_cut
supl_inc$prob_cfor_inc_areaT<-nosf_conf$prob_cfor_inc
supl_inc$are_km2_areaT<-nosf_conf$are_km2
# ammend local tracked colonies
read.csv('C:/seabirds/data/gbr_local_core_foraging_area_inclusion.csv')
supl_inc[supl_inc$mod!='global',]$int<-NA
supl_inc[supl_inc$mod!='global',]$pc<-NA
supl_inc[supl_inc$mod!='global',]$pc_auc<-NA
supl_inc[supl_inc$mod!='global',c(14,16)]<-1 # set all to 100% apart form wtlg heron conf
supl_inc[supl_inc$mod!='global' & supl_inc$md_spgr=='WTLG',16]<-0.87

# summarise and add to tab1
tab1<-cbind(tab1, supl_inc%>%group_by(md_spgr)%>%summarise(mn_incO=min(prob_cfor_inc),me_incO=mean(prob_cfor_inc),mx_incO=max(prob_cfor_inc),
                                         mn_incC=min(prob_cfor_inc_areaT),me_incC=mean(prob_cfor_inc_areaT),mx_incC=max(prob_cfor_inc_areaT)))
#write out tab1
# write.csv(tab1, 'data/conf_area_table_diss2.csv', quote=F, row.names=F)


#write out supl
supl_inc$dsgntn_n<-gsub(',', '@', supl_inc$dsgntn_n)
supl_inc$site_nm<-gsub(',', '@', supl_inc$site_nm)
# write.csv(supl_inc, 'data/gbr_col_corefor_inc_breakdown.csv', quote=F, row.names=F)

# make site level plot main fig (no longer supplement)
corez_nosf$md_spgr<-factor(corez_nosf$md_spgr, levels=c("BRBO", 'MABO', 'RFBO', 'WTST',"WTLG",'FRBD', 'TRBD', "SOTE" , 'NODD', 'TERN'))


p1<-ggplot(data=corez_nosf, aes(x=perc_cut, y=are_km2/1000))+
  geom_segment(data=corez_nosf%>%filter(st_c_ty=='first'), aes(x=0, xend=0.9, y=100, yend=100), colour='purple', linetype='dotted')+
  geom_line(aes(colour=substr(mod, 1, 3),group=site_nm), alpha=0.5)+
  
  geom_segment(data=corez_nosf%>%filter(auc_type=='obs'& auc>0.5), aes(x=perc_cut, xend=perc_cut, y=0, yend=are_km2/1000, group=mod), alpha=0.2, colour='#00b159')+
  geom_segment(data=corez_nosf%>%filter(auc_type=='obs' & auc>0.5), aes(x=0, xend=perc_cut, y=are_km2/1000, yend=are_km2/1000,group=mod), alpha=0.2,colour='#00b159')+
  
  geom_segment(data=corez_nosf%>%filter(st_c_ty=='first'), aes(x=perc_cut, xend=perc_cut, y=0, yend=are_km2/1000, group=mod), alpha=0.2, colour='#d11141')+
  geom_segment(data=corez_nosf%>%filter(st_c_ty=='first'), aes(x=0, xend=perc_cut, y=are_km2/1000, yend=are_km2/1000,group=mod), alpha=0.2, colour='#d11141')+
  
  geom_point(data=filter(corez_nosf, auc==0.5), aes(colour=substr(mod, 1, 3)))+
  geom_point(aes(colour=st_c_ty))+
  geom_point(data=filter(corez_nosf, auc==0.5), aes(colour=substr(mod, 1, 3)), shape=1)+
  facet_wrap(~md_spgr, scales='free', ncol=2)+
  scale_colour_manual(values = c('glo'='#2c56fd','loc'='#ff830f','obs'='#00b159', 'first'='#d11141', 'sim'=NA ))+
  scale_x_continuous(limits=c(0, 0.9), breaks=c(0, 0.25, 0.5, 0.75, 0.9),
                     labels=c('0', '0.25', '0.5', '0.75', '0.9'))+theme_bw()+
  theme(legend.position = "none", panel.grid.minor = element_blank())+
  xlab('Predicted habitat suitability threshold percentile')+
  ylab('Foraging area (thousands of km2)')
#ggsave(p1,  width =4 , height =11.4, units='in',
#       filename='data/modelling/plots/core_cost_confidence_perc_cut_indiv.png')

# Pieces of code to  make mapped part of main figure
## 1)split corez shapefile in many for export for QGIS ##

for(i in unique (corez$md_spgr))
{
  st_write(filter(corez, md_spgr==i & auc==0.5)%>%group_by(dsgntn_n)%>%
             summarise(geometry=st_union(geometry))%>%st_cast(),
           paste0('data/GIS/QGIS_fig_plots/', i, '_rad.shp'), delete_dsn=T)
  st_write(filter(obs_diss, md_spgr==i ),
           paste0('data/GIS/QGIS_fig_plots/', i, '_obs.shp'), delete_dsn=T)
  if('first'%in%corez[corez$md_spgr==i,]$st_c_ty){
    st_write(corez%>%filter(st_c_ty=='first' & md_spgr==i)%>%summarise(geometry=st_union(geometry)),
             paste0('data/GIS/QGIS_fig_plots/', i, '_conf.shp'), delete_dsn=T)}
  print(i)
}

## 2) Mosaic global predictions with local predictions (within obs for rad) for tracked site ##

for_rad<-read_sf('data/GIS/foraging_radii.shp')

for_rad<-for_rad%>%filter(rd_clss=='obs')

for_rad<-for_rad%>%group_by(md_spgr, dsgntn_n)%>%summarise(geometry=st_union(geometry))

pred_list<-list.files('data/modelling/GBR_preds/selected_preds', full.names=T)
mod_pred<-stack(pred_list)

for( i in unique(for_rad$md_spgr))
{
  col<-filter(for_rad, md_spgr==i)
  glob_ras<-subset(mod_pred, paste0(i, '_MultiCol'))
  
  for(j in unique(col$dsgntn_n))
  {
    col_sp<-filter(col, dsgntn_n==j)
    rm(sp_ras1)
    
    if(i=='BRBO' & col_sp$dsgntn_n=='Raine Island, Moulter and MacLennan cays KBA'){
      sp_ras1<-subset(mod_pred, 'BRBO_Raine')}
    if(i=='BRBO' & col_sp$dsgntn_n=='Swain Reefs KBA'){
      sp_ras1<-subset(mod_pred, 'BRBO_Swains')}
    if(i=='MABO' & col_sp$dsgntn_n=='Swain Reefs KBA'){
      sp_ras1<-subset(mod_pred, 'MABO_Swains')}
    if(i=='WTST' & col_sp$dsgntn_n=='Capricornia Cays KBA'){
      sp_ras1<-subset(mod_pred, 'WTST_Heron')}
    if(i=='WTLG' & col_sp$dsgntn_n=='Capricornia Cays KBA'){
      sp_ras1<-subset(mod_pred, 'WTLG_Heron')}
    if(i=='NODD' & col_sp$dsgntn_n=='Capricornia Cays KBA'){
      sp_ras1<-subset(mod_pred, 'NODD_Heron')}
    
    sp_ras2<-crop(sp_ras1, extent(col_sp)) # drop size
    
    sum_rad<-mask(sp_ras2, as(col_sp, 'Spatial')) # mask set to NA
    
    if(which(j==unique(col$dsgntn_n))==1){
      mos_ras <- merge(sum_rad, glob_ras)}else # merge instead of mosaic, first arguement goes over second
      {mos_ras <- merge(sum_rad, mos_ras)}
  }
  print(i)
  
  writeRaster(mos_ras, paste0('C:/seabirds/data/modelling/GBR_preds/glob_local_merge/glob_loc_',i,'.tif'))    
  rm(mos_ras)
}

#### ~~~~ **** ~~~~ ####

#### ~~~~ Make most appropriate refined foraging circle hotspot map for GBR, Fig 6 ~~~~ ####
## Pull optimum refined radii for each sp and combine into 1 shpfile, also calc overap raster

comb_layer<-rbind(
  st_read('data/GIS/QGIS_fig_plots/BRBO_obs.shp'),
  st_read('data/GIS/QGIS_fig_plots/MABO_obs.shp'),
  st_read('data/GIS/QGIS_fig_plots/FRBD_obs.shp'),
  st_read('data/GIS/QGIS_fig_plots/TRBD_obs.shp'),
  st_read('data/GIS/QGIS_fig_plots/TERN_obs.shp'),
  st_read('data/GIS/QGIS_fig_plots/NODD_obs.shp'),
  st_read('data/GIS/QGIS_fig_plots/SOTE_obs.shp'))

rfbo1<-st_read('data/GIS/QGIS_fig_plots/RFBO_conf.shp')
rfbo1$md_spgr<-'RFBO'
rfbo1$are_km2<-999
rfbo1$FID<-NULL

wtst1<-st_read('data/GIS/QGIS_fig_plots/WTST_obs_conf_combined.shp')
wtst1$md_spgr<-'WTST'; wtst1$path<-NULL;wtst1$layer<-NULL;wtst1$FID<-NULL

wtlg1<-st_read('data/GIS/QGIS_fig_plots/WTLG_obs_conf_combined.shp')
          
     
comb_layer<-rbind(comb_layer, rfbo1, wtst1, wtlg1)

#st_write(comb_layer, 'data/GIS/QGIS_fig_plots/all_sp_optimal_refined2.shp')

#count overlap in raster
r1<-raster('data/modelling/GBR_preds/selected_preds/TRBD_MultiCol.tif')
r_sp<-rasterize(as(comb_layer, 'Spatial'), r1, field=1, fun='count')
writeRaster(r_sp, 'data/GIS/QGIS_fig_plots/all_sp_optimal_refined_overlap2.tif') 
            
#### ~~~~ **** ~~~~ ####

####~~~~ return cluster and heatmap validation plot function for Fig 2 ~~~~####
mkVal<-function(my.sp='BRBO', my.metric='AUC', calc.niche=F, text_size=2)
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
  
  #recalc colouring due to colony self change
  my.aucz$auctemp<-my.aucz$auc
  my.aucz[my.aucz$spcol=='MEAN',]$auctemp<-NA
  my.aucz$auc_bin<-cut(my.aucz$auctemp, c(0, 0.6, 0.7, 0.8, 0.9, 1),include.lowest =T)
  my.aucz$auc_bin<-factor(my.aucz$auc_bin, levels = levels(addNA(my.aucz$auc_bin)),
                       labels = c("#cccccc", "#ffffb2", "#fecc5c", '#fd8d3c', '#f03b20', '#c6dbef'), exclude = NULL)
  
  my.aucz$tsstemp<-my.aucz$TSS
  my.aucz[my.aucz$spcol=='MEAN',]$tsstemp<-NA
  my.aucz$tss_bin<-cut(my.aucz$tsstemp, c(0, 0.2, 0.6, 1), include.lowest =T)
  my.aucz$tss_bin<-factor(my.aucz$tss_bin, levels = levels(addNA(my.aucz$tss_bin)),
                       labels = c("#cccccc", "#fecc5c", '#f03b20', '#c6dbef'), exclude = NULL)
  
  
  my.aucz$Resample<-as.character(my.aucz$Resample)
  my.aucz<-filter(my.aucz, !Resample %in% c('EnsembleNrm','EnsembleRaw'))
  my.aucz[my.aucz$Resample=='MultiCol',]$Resample<-'Multi-colony'
  my.aucz$Resample<-factor(my.aucz$Resample) 
  my.aucz$spcol<-factor(my.aucz$spcol)
  my.aucz$Resample<-factor(my.aucz$Resample, levels=c("Multi-colony",paste(hc1$labels[hc1$order])))
  my.aucz$spcol<-factor(my.aucz$spcol,levels=c("MEAN", paste(hc1$labels[hc1$order])))
  
 
  my.aucz_comp<-filter(my.aucz, as.character(Resample)==as.character(spcol) | Resample=='Multi-colony')
  my.aucz_comp<-my.aucz_comp%>%filter(spcol!='MEAN')
  my.aucz_comp$idtemp<-1
  my.aucz_comp[my.aucz_comp$Resample=='Multi-colony',]$idtemp<-2
  my.aucz_comp<-my.aucz_comp%>%arrange(idtemp)%>%group_by(spcol)%>%
    mutate(self_auc=first(auc), mult_auc=last(auc), self_tss=first(TSS), mult_tss=last(TSS))
  
  my.aucz_comp$self_mult_diff_auc<-ifelse(my.aucz_comp$Resample!='Multi-colony',
          my.aucz_comp$auc-my.aucz_comp$mult_auc, my.aucz_comp$auc-my.aucz_comp$self_auc)
  
  my.aucz_comp$self_mult_dclass_auc<-ifelse(abs(my.aucz_comp$self_mult_diff_auc)< 0.05, 'black',
                                        ifelse(my.aucz_comp$self_mult_diff_auc< 0, 'firebrick', 'springgreen4'))
  
  my.aucz_comp$self_mult_diff_tss<-ifelse(my.aucz_comp$Resample!='Multi-colony',
           my.aucz_comp$TSS-my.aucz_comp$mult_tss, my.aucz_comp$TSS-my.aucz_comp$self_tss)
  
  my.aucz_comp$self_mult_dclass_tss<-ifelse(abs(my.aucz_comp$self_mult_diff_tss)< 0.05, 'black',
                                        ifelse(my.aucz_comp$self_mult_diff_tss< 0, 'firebrick', 'springgreen4'))

  mytext<-filter(my.aucz, as.character(Resample)!=as.character(spcol) & Resample!='Multi-colony')
  mytext<-rbind(mytext,filter(my.aucz, Resample=='Multi-colony'&spcol=='MEAN'))
  
  if(my.metric=='AUC')
  {  
  val_plot<-ggplot(my.aucz, aes(x = spcol, y = Resample)) + 
    geom_raster(aes(fill=auc_bin)) +scale_fill_identity()+ 
    geom_text(data=mytext, aes(label=round(auc, 2)*100), size=text_size)+
    geom_text(data=my.aucz_comp, aes(label=round(auc, 2)*100,colour=self_mult_dclass_auc), size=text_size, fontface='bold')+
    scale_colour_identity()+theme_bw() + 
    theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))+ylab('Predictions from')+xlab('Predicting to')}
  
  if(my.metric=='TSS')
  {
    val_plot<-ggplot(my.aucz, aes(x = spcol, y = Resample)) + 
      geom_raster(aes(fill=tss_bin)) +scale_fill_identity()+  
      geom_text(data=mytext,aes(label=round(TSS, 2)*100), size=text_size)+
      geom_text(data=my.aucz_comp, aes(label=round(TSS, 2)*100,colour=self_mult_dclass_tss), size=text_size, fontface='bold')+
      scale_colour_identity()+theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.3),
                         axis.text.y=element_text(size=9),
                         plot.title=element_text(size=11))+ylab('Predictions from')+xlab('Predicting to')}
  
  if(calc.niche==T){
  # niche overlap
  dat<-read.csv(paste0('data/modelling/kernhull_pts_sample/', my.sp, '_kernhull_sample.csv'))
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

#### ~~~~ Make seperate figure componenets for Fig 2, for export to CorelDraw ~~~~####
#https://cran.r-project.org/web/packages/ggplotify/vignettes/ggplotify.html

# BRBO 
brbo_auc<-mkVal('BRBO', 'AUC', calc.niche = F, text_size = 2.5)
brbo_tss<-mkVal('BRBO', 'TSS', calc.niche = F, text_size = 2.5)
brbo_auc_dend<-ggdendrogram(data = as.dendrogram(brbo_auc[[1]]), rotate = T)
brbo_tss_dend<-ggdendrogram(data = as.dendrogram(brbo_tss[[1]]), rotate = T)
# MABO
mabo_auc<-mkVal('MABO', 'AUC', calc.niche = F, text_size = 2.5)
mabo_tss<-mkVal('MABO', 'TSS', calc.niche = F, text_size = 2.5)
mabo_auc_dend<-ggdendrogram(data = as.dendrogram(mabo_auc[[1]]), rotate = T)
mabo_tss_dend<-ggdendrogram(data = as.dendrogram(mabo_tss[[1]]), rotate = T)
# RFBO
rfbo_auc<-mkVal('RFBO', 'AUC', calc.niche = F, text_size = 2.5)
rfbo_tss<-mkVal('RFBO', 'TSS', calc.niche = F, text_size = 2.5)
rfbo_auc_dend<-ggdendrogram(data = as.dendrogram(rfbo_auc[[1]]), rotate = T)
rfbo_tss_dend<-ggdendrogram(data = as.dendrogram(rfbo_tss[[1]]), rotate = T)
# FRBD
frbd_auc<-mkVal('FRBD', 'AUC', calc.niche = F, text_size = 2.5)
frbd_tss<-mkVal('FRBD', 'TSS', calc.niche = F, text_size = 2.5)
frbd_auc_dend<-ggdendrogram(data = as.dendrogram(frbd_auc[[1]]), rotate = T)
frbd_tss_dend<-ggdendrogram(data = as.dendrogram(frbd_tss[[1]]), rotate = T)
# TRBD
trbd_auc<-mkVal('TRBD', 'AUC', calc.niche = F, text_size = 2.5)
trbd_tss<-mkVal('TRBD', 'TSS', calc.niche = F, text_size = 2.5)
trbd_auc_dend<-ggdendrogram(data = as.dendrogram(trbd_auc[[1]]), rotate = T)
trbd_tss_dend<-ggdendrogram(data = as.dendrogram(trbd_tss[[1]]), rotate = T)
# WTST
wtst_auc<-mkVal('WTST', 'AUC', calc.niche = F, text_size = 2.5)
wtst_tss<-mkVal('WTST', 'TSS', calc.niche = F, text_size = 2.5)
wtst_auc_dend<-ggdendrogram(data = as.dendrogram(wtst_auc[[1]]), rotate = T)
wtst_tss_dend<-ggdendrogram(data = as.dendrogram(wtst_tss[[1]]), rotate = T)
# WTLG
wtlg_auc<-mkVal('WTLG', 'AUC', calc.niche = F, text_size = 2.5)
wtlg_tss<-mkVal('WTLG', 'TSS', calc.niche = F, text_size = 2.5)
wtlg_auc_dend<-ggdendrogram(data = as.dendrogram(wtlg_auc[[1]]), rotate = T)
wtlg_tss_dend<-ggdendrogram(data = as.dendrogram(wtlg_tss[[1]]), rotate = T)
# NODD
nodd_auc<-mkVal('NODD', 'AUC', calc.niche = F, text_size = 2.5)
nodd_tss<-mkVal('NODD', 'TSS', calc.niche = F, text_size = 2.5)
nodd_auc_dend<-ggdendrogram(data = as.dendrogram(nodd_auc[[1]]), rotate = T)
nodd_tss_dend<-ggdendrogram(data = as.dendrogram(nodd_tss[[1]]), rotate = T)
# SOTE
sote_auc<-mkVal('SOTE', 'AUC', calc.niche = F, text_size = 2.5)
sote_tss<-mkVal('SOTE', 'TSS', calc.niche = F, text_size = 2.5)
sote_auc_dend<-ggdendrogram(data = as.dendrogram(sote_auc[[1]]), rotate = T)
sote_tss_dend<-ggdendrogram(data = as.dendrogram(sote_tss[[1]]), rotate = T)
# TERN
tern_auc<-mkVal('TERN', 'AUC', calc.niche = F, text_size = 2.5)
tern_tss<-mkVal('TERN', 'TSS', calc.niche = F, text_size = 2.5)
tern_auc_dend<-ggdendrogram(data = as.dendrogram(tern_auc[[1]]), rotate = T)
tern_tss_dend<-ggdendrogram(data = as.dendrogram(tern_tss[[1]]), rotate = T)

# AUC plots
brbo<-(brbo_auc[[2]]+ggtitle('A) Brown Booby'))+brbo_auc_dend+
  plot_layout(ncol=2, nrow=1, widths=c(3,1))
ggsave(plot=brbo, filename='C:/seabirds/data/modelling/plots/BRBO.png',width = 7, height =5)

maborfbo<-(mabo_auc[[2]]+ggtitle('A) Masked Booby'))+mabo_auc_dend+(rfbo_auc[[2]]+ggtitle('B) Red-footed Booby'))+rfbo_auc_dend+
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
ggsave(plot=noddsotetern, filename='C:/s eabirds/data/modelling/plots/NODD_SOTE_TERN.png',width = 7, height =12)

# eps plot only
ggsave(plot=brbo_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
       axis.title.x=element_blank(), axis.ticks=element_blank(),axis.title.y=element_blank()), width = 17*0.5, height =17*0.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_BRBO.eps')
ggsave(plot=mabo_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(),  axis.ticks=element_blank(),axis.title.y=element_blank()), width = 13*0.5, height =13*0.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_MABO.eps')
ggsave(plot=rfbo_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(),  axis.ticks=element_blank(),axis.title.y=element_blank()), width = 12*0.5, height =12*0.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_RFBO.eps')
ggsave(plot=frbd_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()), width = 10*0.5, height =10*0.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_FRBD.eps')
ggsave(plot=trbd_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()), width = 12*0.5, height =12*0.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_TRBD.eps')
ggsave(plot=wtst_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(),  axis.ticks=element_blank(),axis.title.y=element_blank()), width = 9*0.5, height =9*0.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_WTST.eps')
ggsave(plot=wtlg_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(), axis.ticks=element_blank(),axis.title.y=element_blank()), width = 8*0.5, height =8*0.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_WTLG.eps')
ggsave(plot=nodd_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(), axis.ticks=element_blank(),axis.title.y=element_blank()), width = 3.6, height =3.6, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_NODD.eps')
ggsave(plot=sote_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(), axis.ticks=element_blank(),axis.title.y=element_blank()), width = 2.5, height =2.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_SOTE.eps')
ggsave(plot=tern_auc[[2]]+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.title.x=element_blank(), axis.ticks=element_blank(),axis.title.y=element_blank()), width = 2.5, height =2.5, units='cm',
       filename='C:/seabirds/data/modelling/plots/matr_TERN.eps')


#### ~~~~ *** ~~~~ ####

####~~~~ local adaptation mantel and plots function, for suppl mat ~~~~####
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

