# Plotting and mapping
library(ggplot2)
#library(ggdendro)
library(dplyr)
#library(grid)
#library(gridExtra)
library(tidyr)
library(ggplotify)
library(patchwork)
library(vegan)
library(sf)

####~~~~ read in data ~~~~####
aucz_out<-read.csv('C:/seabirds/data/mod_validation_vals.csv')
matx_out<-read.csv('C:/seabirds/data/mod_clustering_vals.csv')
t_qual<-read.csv('C:/seabirds/data/tracking_trip_decisions.csv')
colz<-st_read('C:/seabirds/data/GIS/trackingID_colony_locs.shp')


#### sp-col summary table

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
t_tripmetric<-t_qual_ret%>%group_by(sp_group, sp, coly)%>%
  summarise(max_for=max(max_dist), med_for=median(max_dist), breedstage=paste(unique(breedstage), collapse=' '))%>%as.data.frame()
spcol_tab<-data.frame(spcol_tab, t_tripmetric[,c(4,5,6)])
# tidy
spcol_tab[spcol_tab$coly=='chick',]$coly<-'Rat' # edit to name
spcol_tab[spcol_tab$sp_group=='WTST' | spcol_tab$sp_group=='WTLG',]$sp<-'WTSH'

#write.csv(spcol_tab, 'C:/seabirds/data/sp_col_summary.csv', quote=F, row.names=F)

#### sp-group summary table



aucz_out%>%filter(as.character(spcol)!=as.character(Resample) & spcol!='SUM')%>%group_by(sp, Resample)%>%
  summarise(mean_auc=mean(auc), sd_auc=sd(auc))





# AUC/TSS density plot

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

  my.aucz$Resample<-factor(my.aucz$Resample, levels=c("MultiCol", paste(hc1$labels[hc1$order])))
  my.aucz$spcol<-factor(my.aucz$spcol,levels=c("SUM", paste(hc1$labels[hc1$order])))
  
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