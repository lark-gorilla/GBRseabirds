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
# tidy
spcol_tab[spcol_tab$coly=='chick',]$coly<-'Rat' # edit to name
spcol_tab[spcol_tab$sp_group=='WTST' | spcol_tab$sp_group=='WTLG',]$sp<-'WTSH'
spcol_tab[spcol_tab$coly=='Adele' & spcol_tab$sp=='BRBO',]$max_for<-139.6  # manual edit

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

for_rang<-spcol_tab%>%group_by(Species.Group)%>%summarise(min.for=min(Max.roraging.range),
                                                med.for=median(Max.roraging.range),
                                                max.for=max(Max.roraging.range))

sp_col_summr<-aucz_out%>%filter(as.character(spcol)!=as.character(Resample) & spcol!='SUM')%>%group_by(sp, Resample)%>%
  summarise(mean_auc=mean(auc), sd_auc=sd(auc), max_auc=max(auc),
            mean_tss=mean(TSS), sd_tss=sd(TSS),max_tss=max(TSS))

sp_col_summr<-sp_col_summr%>%ungroup()%>%
  group_by(sp)%>%mutate(auc_rank=rank(-mean_auc, 'first'),tss_rank=rank(-mean_tss, 'first'))

mean_valz<-sp_col_summr%>%ungroup()%>%filter(Resample!='MultiCol')%>%group_by(sp)%>%
  summarise(mn_auc=mean(mean_auc), sd_auc=sd(mean_auc),
            mn_max_auc=mean(max_auc), sd_max_auc=sd(max_auc),
            mn_tss=mean(mean_tss), sd_tss=sd(mean_tss),
            mn_max_tss=mean(max_tss), sd_max_tss=sd(max_tss))

multi_rank<-sp_col_summr%>%filter(Resample=='MultiCol')%>%select(auc_rank, tss_rank)

#GBR pred
aucz_out%>%filter(as.character(spcol)!=as.character(Resample) & spcol%in%c('Raine', 'Swains', 'Heron')& auc>0.69)%>%group_by(sp, Resample)
GBR_pred<-sp_col_summr%>%ungroup()%>%filter(Resample!='MultiCol')%>%group_by(sp)

bind_out<-bind_cols(mean_valz, multi_rank, for_rang)

bind_out$mn_sumr<-paste0(round(bind_out$mn_auc, 2),'±', round(bind_out$sd_auc, 2))
bind_out$mx_sumr<-paste0(round(bind_out$mn_max_auc, 2),'±', round(bind_out$sd_max_auc, 2))
bind_out$mn_sumr_ts<-paste0(round(bind_out$mn_tss, 2),'±', round(bind_out$sd_tss, 2))
bind_out$mx_sumr_ts<-paste0(round(bind_out$mn_max_tss, 2),'±', round(bind_out$sd_max_tss, 2))


bind_out2<-bind_out%>%select(sp, mn_sumr, mx_sumr, mn_sumr_ts, mx_sumr_ts, auc_rank, tss_rank, min.for, med.for, max.for)
#write.csv(bind_out2, 'C:/seabirds/data/sp_main_summary.csv', quote=F, row.names=F)

#### ~~~~ **** ~~~~ ####

#### ~~~~ Make GBR colony radii ~~~~####
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

#### ~~~~ Make plots ~~~~ ####
#read in spatial data
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

  p1<-ggplot() +
  layer_spatial(data=r1) +
  geom_sf(data=filter(colz, md_spgr==spg), colour='yellow', size=0.5) +
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
                     option='magma', na.value = NA)+
  theme(legend.position = 'none')
  return(p1)}

#theme(legend.position = c(0.2, 0.3), legend.background = element_rect(fill = "grey"),
 #     legend.key = element_rect(fill = "grey"))
# check higher res/different png dimensions
# patchwork <-p1 + p2 + p3 + p4 + plot_layout(ncol = 4,guides = "collect")

#### ~~~~ **** ~~~~ #####
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
(p_brbo+ggtitle('A) Brown Booby')+p_mabo+ggtitle('B) Masked Booby')+
 p_rfbo+ggtitle('C) Red-footed Booby')+p_frbd+ggtitle('D) Frigatebird species-group')+
  plot_layout(ncol=2, nrow=2))
dev.off()

png(paste0('C:/seabirds/outputs/maps/gbr_wide/tropbd2wtsh2sote.png'),width = 8.3, height =11.7 , units ="in", res =600)
(p_trbd+ggtitle('E) Tropicbird species-group')+p_wtst+ggtitle('F) Wedge-tailed Shearwater short trips')+
    p_wtlg+ggtitle('G) Wedge-tailed Shearwater long trips')+p_sote+ggtitle('H) Sooty Tern')+
    plot_layout(ncol=2, nrow=2))
dev.off()

png(paste0('C:/seabirds/outputs/maps/gbr_wide/nodd2tern.png'),width = 8.3, height =5.85 , units ="in", res =600)
(p_nodd+ggtitle('I) Noddy species-group')+p_tern+ggtitle('J) Tern species-group')+
    plot_layout(ncol=2))
dev.off()

#### ~~~~ kba/site local plot function ~~~~ ####
mk_kbaplot<-function(site="Capricornia Cays KBA"){
  
  rad_diss_site<-filter(for_rad, dsgntn_n==site)%>%group_by(md_spgr, rd_clss)%>%
    summarize(spcol=first(spcol),geometry = st_union(geometry))
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
        breaks = unique(rad_diss_site$spcol), guide = "legend")+ggtitle(paste0('A) ', site))
  
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
    layer_spatial(data=r1) +
    geom_sf(data=filter(gbr_reef, FEAT_NAME=='Reef'), fill='NA', colour=alpha('black',0.5))+
    geom_sf(data=filter(colz, dsgntn_n==site& md_spgr==i), shape = 23, fill = "darkred") +
    geom_sf(data=filter(rad_diss_site, md_spgr==i), aes(colour=rd_clss), fill='NA') +
    geom_sf(data=gbrmp, col='white', fill='NA') +
    geom_sf(data=land, col='black', fill='grey') +
    theme_bw()+
    annotation_scale(location = "bl")+  
    annotation_north_arrow(location = "tr", which_north = "true")+
    labs(x='Longitude', y='Latitude')+
    coord_sf(xlim = plim[c(1,3)], ylim = plim[c(2,4)], expand = T)+
    scale_colour_manual('Forgaing\nradii', values=c('#00FFFF','#66FFCC', '#00FF66' ), labels=c(
      'Max', 'Med', 'Min'))+
    scale_fill_viridis('Likely\nforaging\nhabitat',option='magma',
                       breaks=c(mn, mx), labels=c('low', 'high'),na.value = NA)+
    ggtitle(paste0(LETTERS[(which(i== unique(rad_diss_site$md_spgr)))+1], ') ', i))
  
  if(which(i== unique(rad_diss_site$md_spgr))!=1){p2<-p2+theme(legend.position = 'none')}
  p_spz[[which(i== unique(rad_diss_site$md_spgr))]]<-p2
  print(i)
  }
  
  pw_plot<-p1+p_spz+ plot_layout(ncol=3, guides='collect')
  
  return(pw_plot)}
#### ~~~~ **** ~~~~ #####

#### ~~~~ Make local KBA  plots ~~~~ ####

for(k in unique(for_rad$dsgntn_n))
{
  multip<-mk_kbaplot(site=k) 
  png(paste0('C:/seabirds/outputs/maps/local_site/',gsub( ',', '',k),'.png'),width = 8.3, height =11.7 , units ="in", res =600)
  multip
  dev.off()
}

#### ~~~~ **** ~~~~ #####

####~~~~ AUC/TSS density plot ~~~~####

densdat<-aucz_out%>%filter(as.character(spcol)!=as.character(Resample) & spcol!='SUM')%>%
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
mkVal<-function(my.sp='BRBO', my.metric='AUC', calc.niche=T)
{
  #https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html
  my.matx<-matx_out[matx_out$sp==my.sp,]
  my.aucz<-aucz_out[aucz_out$sp==my.sp,]
  
  if(my.metric=='AUC')
  {  
  d1<-with(my.matx[my.matx$Resample!=my.matx$spcol,c(10,11,5)], 
           structure(auc, Size = length(unique(my.matx$Resample)),
                     Labels = unique(my.matx$Resample),
                     Diag = F, Upper = FALSE,method = "user", class = "dist"))}
  if(my.metric=='TSS')
    {
    d1<-with(my.matx[my.matx$Resample!=my.matx$spcol,c(10,11,9)], 
             structure(TSS, Size = length(unique(my.matx$Resample)),
                       Labels = unique(my.matx$Resample),
                       Diag = F, Upper = FALSE,method = "user", class = "dist"))}
  
  d1<-(2-d1)
  hc1<-hclust(d1, method='average')
  hc_dend<-ggdendrogram(data = as.dendrogram(hc1), rotate = T)
  
  # make sum/ave vals
  my.aucz<-my.aucz[,c(1:10, 16, 18)]
  my.aucz<-filter(my.aucz, spcol!='SUM')
  my.aucz$spcol<-factor(my.aucz$spcol)
  
  temp1<-my.aucz%>%group_by(Resample)%>%
    filter(as.character(spcol)!=as.character(Resample))%>%summarise_if(is.numeric ,mean)
  my.aucz<-rbind(my.aucz,data.frame(sp=my.sp,temp1[,1], spcol='MEAN', temp1[,2:8], auc_bin='#c6dbef', tss_bin='#c6dbef'))

  my.aucz$Resample<-factor(my.aucz$Resample, levels=c("MultiCol", paste(hc1$labels[hc1$order])))
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
sp_auc<-mkVal('BRBO', 'AUC', calc.niche = T)
sp_tss<-mkVal('BRBO', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/BRBO_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
    plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# MABO
sp_auc<-mkVal('MABO', 'AUC', calc.niche = T)
sp_tss<-mkVal('MABO', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/MABO_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# RFBO
sp_auc<-mkVal('RFBO', 'AUC', calc.niche = T)
sp_tss<-mkVal('RFBO', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/RFBO_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# FRBD
sp_auc<-mkVal('FRBD', 'AUC', calc.niche = T)
sp_tss<-mkVal('FRBD', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/FRBD_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# TRBD
sp_auc<-mkVal('TRBD', 'AUC', calc.niche = T)
sp_tss<-mkVal('TRBD', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/TRBD_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# WTST
sp_auc<-mkVal('WTST', 'AUC', calc.niche = T)
sp_tss<-mkVal('WTST', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/WTST_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# WTLG
sp_auc<-mkVal('WTLG', 'AUC', calc.niche = T)
sp_tss<-mkVal('WTLG', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/WTLG_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# NODD
sp_auc<-mkVal('NODD', 'AUC', calc.niche = T)
sp_tss<-mkVal('NODD', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/NODD_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# SOTE
sp_auc<-mkVal('SOTE', 'AUC', calc.niche = T)
sp_tss<-mkVal('SOTE', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/SOTE_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()

# TERN
sp_auc<-mkVal('TERN', 'AUC', calc.niche = T)
sp_tss<-mkVal('TERN', 'TSS', calc.niche = F)

p1<-ggdendrogram(data = as.dendrogram(sp_auc[[1]]), rotate = T)
p2<-ggdendrogram(data = as.dendrogram(sp_tss[[1]]), rotate = T)

png(paste0('C:/seabirds/data/modelling/plots/TERN_val.png'),width = 7, height =10 , units ="in", res =300)
(sp_auc[[2]]+ggtitle('A)'))+p1+(sp_tss[[2]]+ggtitle('B)'))+p2+
  plot_layout(ncol=2, nrow=2, widths=c(3,1))
dev.off()





png(paste0('C:/seabirds/data/modelling/plots/temp.png'),width = 8, height =6 , units ="in", res =600)
grid.newpage()
print(brbo_auc[[2]], vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(brbo_auc[[1]], vp = viewport(x = 0.90, y = 0.58, width = 0.2, height = 0.82))
dev.off()