# Modelling (kernelUD approach)

library(dplyr)
library(ggplot2)
#library(GGally)
library(ranger)
library(caret)
library(pROC)
library(vegan)
library(raster)
library(ncf)
library(gridExtra)

set.seed(24) # so that random samples and random forest models can be reproduced
normalized<-function(x){(x-min(x))/(max(x)-min(x))} # normalise (0-1) function

####~~ Random sampling of psuedo absences fixed for reproducability amd data compression ~~~####
#sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG','FRBD', 'TRBD', 'NODD', 'TERN')
#for(k in sp_groups)
#{
#  dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts/', k, '_kernhull.csv'))
#  names(dat)[names(dat)=='layer']<-'forbin'
#  dat$forbin<-as.factor(dat$forbin)
#  levels(dat$forbin)<-c("PsuedoA", "Core")
#  dat[dat$forbin=='Core',]$weight=0 # change from NA for no.omit
  
  # set bathy vals >0 to NA then make positive
#  dat[dat$bth>0,]$bth<-NA
#  dat$bth<-sqrt(dat$bth^2)
#  dat<-na.omit(dat) # so pA's aren't assigned to land

  # RM sp name unless k has conflicts not FRBD or TERN
#  if(k %in% c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
#              'TRBD', 'NODD')){dat$spcol<-substr(dat$spcol, 6, nchar(as.character(dat$spcol)))}
#  dat$spcol<-factor(dat$spcol)
  #dat%>%group_by(spcol)%>%summarise(nP=length(which(forbin=='Core')),nA=length(which(forbin=="PsuedoA"))) # don't need
#    dat<-dat%>%group_by(spcol)%>%mutate(npres=length(which(forbin=='Core')))

    # 3:1 Pa to A, sampled wothout replacement, over normalised inv coldist surface
#  dat<-rbind(dat%>%filter(forbin=='Core'), dat%>%filter(forbin=="PsuedoA")%>%group_by(spcol)%>%
#               mutate(w2=normalized(weight))%>% 
#               sample_n(size=(unique(npres)*3), replace=F, weight=w2))%>%as.data.frame() 
#  dat$npres<-NULL
#  dat$w2<-NULL
#  write.csv(dat,paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', k, '_kernhull_sample.csv'), row.names=F, quote=F)
#  print(k)
#}
####~~~ * ~~~####

####~~~ Loop above creates data that are handed to modelling_kernUD.R for tuning on HPC ~~~####

####~~~ Visualise RF tuning results to manually select optimum hyperparameters ~~~####
#lsf<-list.files('C:/seabirds/data/modelling/rf_tuning')
#indiv_tune<-NULL; for(i in lsf[grep('indiv',lsf)])
#{indiv_tune<-rbind(indiv_tune, data.frame(read.csv(
#  paste0('C:/seabirds/data/modelling/rf_tuning/', i)), sp=substr(i, 1,4)))}
#all_tune<-NULL; for(i in lsf[grep('all',lsf)])
#{all_tune<-rbind(all_tune, data.frame(read.csv(
#  paste0('C:/seabirds/data/modelling/rf_tuning/', i)), sp=substr(i, 1,4)))}
#sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG', 'FRBD', 'TRBD', 'NODD', 'TERN')
#for(i in sp_groups)
#{
#  p1<-ggplot(data=filter(all_tune, sp==i), aes(x=mtry, colour=factor(min.node.size)))+
#    geom_line(aes(y=auc))+geom_line(aes(y=TSS),linetype='dashed')+
#    geom_point(aes(y=auc), shape=1)+geom_point(aes(y=TSS), shape=2)+facet_wrap(~Resample)+theme(legend.position = "none")  
#  png(paste0('C:/seabirds/data/modelling/rf_tuning/plots/',i, '_all_col.png'),width = 6, height =6 , units ="in", res =600)
#  print(p1)
#  dev.off()
#  for(j in unique(filter(indiv_tune, sp==i)$Resample)){
#    p2<-ggplot(data=filter(indiv_tune, sp==i & Resample==j), aes(x=mtry, colour=factor(min.node.size)))+
#    geom_line(aes(y=auc))+geom_line(aes(y=TSS),linetype='dashed')+
#    geom_point(aes(y=auc), shape=1)+geom_point(aes(y=TSS), shape=2)+facet_wrap(~spcol)+theme(legend.position = "none")  
#    png(paste0('C:/seabirds/data/modelling/rf_tuning/plots/',i, 'indiv_', j,'.png'),width = 6, height =6 , units ="in", res =600)
#    print(p2)
#    dev.off()}
#  print(i)
#}
####~~~ * ~~~####

# Read in optimal hyperparameters
my_hyp<-read.csv('C:/seabirds/data/rf_optimal_hyp.csv')

#read in GBR ocean data
gbr<-read.csv('C:/seabirds/data/pred_area_modelready_2km.csv')
#rename
names(gbr)<-c('x', 'y', 'sst', 'sst_sd', 'chl', 'chl_sd',
              'mfr', 'mfr_sd', 'pfr', 'pfr_sd', 'bth', 'slp')
# and clean
gbr[gbr$bth>0,]$bth<-0
gbr$bth<-sqrt(gbr$bth^2)
# Remember nearshore front values which ==0 should be NA


# read in 2km rasterize template
templ<-raster('C:/seabirds/data/GIS/pred_area_ras_template2km.tif')

# read in data
sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

matx_out<-NULL # capture hclust for plotting
aucz_out<-NULL # capture validation for plotting

for(k in sp_groups)
{
  
  dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', k, '_kernhull_sample.csv'))
  dat$X<-NULL

  # Read in tuning results
  indiv_col_tune<-read.csv(paste0('C:/seabirds/data/modelling/rf_tuning/',k, '_indiv_col_tune.csv'))
  all_col_tune<-read.csv(paste0('C:/seabirds/data/modelling/rf_tuning/',k, '_all_col_tune.csv'))

  # lookup optimal vals
  indiv_col_tune<-na.omit(left_join(indiv_col_tune, filter(my_hyp, sp==k),
                                    by=c('Resample', 'mtry', 'min.node.size')))
  indiv_col_tune$sp<-NULL
  all_col_tune<-na.omit(left_join(all_col_tune, filter(my_hyp, Resample=='MultiCol' & sp==k),
                            by=c('mtry', 'min.node.size')))
  
  if(k=='RFBO')
  {
    indiv_col_tune<-filter(indiv_col_tune, Resample!='Christmas') 
    indiv_col_tune<-filter(indiv_col_tune, spcol!='Christmas') 
    all_col_tune<-filter(all_col_tune, Resample.x!='Christmas') 
    dat<-filter(dat, spcol!='Christmas') 
  }


sp_store<-NULL
var_imp<-NULL
for( i in unique(dat$spcol))
{
  #pr1<-table(dat[dat$spcol==i,]$forbin)[1]/table(dat[dat$ID==i,]$forbin)[2]

  rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp,
                data=dat[dat$spcol==i,], num.trees=500, 
                mtry=unique(filter(indiv_col_tune, Resample==i)$mtry), 
                min.node.size=unique(filter(indiv_col_tune, Resample==i)$min.node.size),
                splitrule = "gini",  importance='impurity',probability =T, seed=24)  
    
  # predict to other colonies 
  dat$p1<-predict(rf1, data=dat)$predictions[,1] # prob of foraging 0-1
  names(dat)[which(names(dat)=='p1')]<-paste(k, i, sep='_')
  
  # predict to training and add to validation data.frame
  p1<-predict(rf1, data=dat[dat$spcol==i,])$predictions[,1] # prob of foraging 0-1 Note column 1 == Core
  auc1<-pROC::roc(dat[dat$spcol==i,]$forbin, p1, levels=c('PsuedoA', 'Core'), direction="<")
  co1<-coords(pROC::roc(dat[dat$spcol==i,]$forbin, p1, levels=c('PsuedoA', 'Core'), direction="<"),'best', best.method='youden', transpose=F)
  
  indiv_col_tune<-rbind(indiv_col_tune, 
  data.frame(Resample=i, spcol=i,  mtry=rf1$mtry, min.node.size=rf1$min.node.size, auc=as.double(auc1$auc), 
             thresh=co1$threshold[1], sens=co1$sensitivity[1], spec=co1$specificity[1],
             TSS=co1$sensitivity[1]+co1$specificity[1]-1))

 
  #varib importance
  var_imp<-rbind(var_imp, data.frame(spcol=i, t(ranger::importance(rf1)/max(ranger::importance(rf1)))))
  
  # predict to GBR
  gbr$p1<-predict(rf1, data=gbr)$predictions[,1] # prob of foraging 0-1
  names(gbr)[which(names(gbr)=='p1')]<-paste(k, i, sep='_')
  
  # predict to GBR training data for certain sp
  #gbr_valdat$p1<-predict(rf1, data=gbr_valdat)$predictions[,1]
  #names(gbr_valdat)[which(names(gbr_valdat)=='p1')]<-paste(k, i, sep='_')
  
  # Internval cross validation for GBR models
  #coltemp<-dat[dat$spcol==i,]
  #coltempcore<-coltemp[coltemp$forbin=='Core',]
  #coltempcore$clust<-kmeans(coltempcore[,14:15], 5)$cluster
  
  #coltemp_subs<-NULL
  #for(h in 1:max(coltempcore$clust))
  #{
  #maxn<-nrow(coltemp[coltemp$forbin=='PsuedoA',])
  #PA_samp<-coltemp[coltemp$forbin=='PsuedoA',][sample(1:maxn, (nrow(coltempcore[coltempcore$clust==h,])*3), replace=F),]
  #PA_samp$clust<-h
  #coltemp_subs<-rbind(coltemp_subs, rbind(coltempcore[coltempcore$clust==h,], PA_samp))
  #}
# output sample so it is fixed 

# subs_preds<-NULL
#  for(m in 1:max(coltempcore$clust))
#  {
#  rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp,
#              data=coltemp_subs[coltemp_subs$clust==m,], num.trees=500, 
#              mtry=unique(filter(indiv_col_tune, Resample==i)$mtry), 
#              min.node.size=unique(filter(indiv_col_tune, Resample==i)$min.node.size),
#              splitrule = "gini",  importance='impurity',probability =T, seed=24)
#  
#  subs_preds<-cbind(subs_preds, predict(rf1, data=dat)$predictions[,1])
#  #print(head(subs_preds))
#  print(m)
#  }

# predict to other colonies 
#dat$p1<-rowSums(subs_preds) # prob of foraging 0-1
#names(dat)[which(names(dat)=='p1')]<-paste(k, i, sep='_')  

  # messing with stratified sampling
#model2 = randomForest(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp, 
#                                           data = coltemp_subs, 
#                                          sampsize=rep(min(table(coltemp_subs$clust)),5),
#                      strata=coltemp_subs$clust)
#dat$p1<-predict(model2, newdata=dat, type='prob')[,1]
#names(dat)[which(names(dat)=='p1')]<-paste(k, i, sep='_')  

  
  #SPAC assessment
  #wtdist<-max(dat[dat$spcol==i,]$weight)-100 # 100 km from col
  #residz<-as.integer(as.character(dat[dat$spcol==i & dat$weight>wtdist,]$forbin))-
  #  dat[dat$spcol==i & dat$weight>wtdist,i]
  #sp1<-spline.correlog(dat[dat$spcol==i& dat$weight>wtdist,]$Longitude,
  #                     dat[dat$spcol==i& dat$weight>wtdist,]$Latitude, residz,
  #                     na.rm=T,latlon=T,resamp=10)
  #plot(sp1)
  #sp_store<-rbind(sp_store, data.frame(spcol=i, Dist=sp1$boot$boot.summary$predicted$x[1,],
  #                     SPAC_025=sp1$boot$boot.summary$predicted$y[3,],
  #                     SPAC_Ave=sp1$boot$boot.summary$predicted$y[6,],
  #                     SPAC_95=sp1$boot$boot.summary$predicted$y[9,]))
  print(i)
  print(Sys.time()) 
  rm(rf1) # save space
}

# make full model and predict to GBR
# check if adding case weights improves(https://github.com/imbs-hl/ranger/issues/167)
rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp ,
            data=dat, num.trees=500,  
            mtry=unique(all_col_tune$mtry), 
            min.node.size=unique(all_col_tune$min.node.size),
            splitrule = "gini",  importance='impurity',probability =T, seed=24)

gbr$MultiCol<-predict(rf1, data=gbr)$predictions[,1]
names(gbr)[which(names(gbr)=='MultiCol')]<-paste(k, 'MultiCol', sep='_')

# predict to GBR training data for certain sp
#rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp ,
#            data=dat[dat$spcol!='Heron',], num.trees=500,  
#            mtry=unique(all_col_tune$mtry), 
#            min.node.size=unique(all_col_tune$min.node.size),
#            splitrule = "gini",  importance='impurity',probability =T, seed=24)
#gbr_valdat$p1<-predict(rf1, data=gbr_valdat)$predictions[,1]
#names(gbr_valdat)[which(names(gbr_valdat)=='p1')]<-paste(k, 'MultiCol_Heron', sep='_')
# Write out predictions for sp with GBR-local data for validation plots
#write.csv(gbr_valdat, 'C:/seabirds/data/modelling/GBR_validation/NODD_gbr_val.csv', quote=F, row.names = F)


# Write out spatial predictions
spdf<-SpatialPointsDataFrame(SpatialPoints(gbr[,1:2], proj4string = CRS(projection(templ))),
                             data=gbr[,grep(k, names(gbr))])

for(j in 1:length(spdf@data))
{
  p1<-rasterize(spdf, templ, field=names(spdf@data)[j])
  writeRaster(p1, paste0('C:/seabirds/data/modelling/GBR_preds/', gsub(' ', '_', names(spdf@data)[j]), '.tif'), overwrite=T)
  #if(i==1){ps<-p1}else{ps<-stack(ps, p1)}
}

#write.csv(sp_store, paste0('C:/seabirds/data/modelling/SPAC/',k,'_spac.csv'), quote=F, row.names = F)
#write.csv(var_imp, paste0('C:/seabirds/data/modelling/var_imp/',k,'_var_imp.csv'), quote=F, row.names = F)

# Clustering sites based on predictive ability
if(k=='SOTE'){
  indiv_col_tune$Resample<-as.character(indiv_col_tune$Resample)
  indiv_col_tune$spcol<-as.character(indiv_col_tune$spcol)
  indiv_col_tune[indiv_col_tune$Resample=='chick',]$Resample<-'Rat'
  indiv_col_tune[indiv_col_tune$spcol=='chick',]$spcol<-'Rat'
  all_col_tune$Resample.x<-as.character(all_col_tune$Resample.x)
  all_col_tune[all_col_tune$Resample.x=='chick',]$Resample.x<-'Rat'
  dat$spcol<-as.character(dat$spcol)
  dat[dat$spcol=='chick',]$spcol<-'Rat'}

indiv_col_tune$id<-apply(as.data.frame(indiv_col_tune), 1, FUN=function(x){paste(sort(c(as.character(x[1]), as.character(x[2]))), sep='.', collapse='.')})
matxdat<-indiv_col_tune%>%group_by(id)%>%summarise_if(is.numeric, sum)
indiv_col_tune$id<-NULL
matxdat$Resample<-unlist(lapply(strsplit(matxdat$id, '\\.'), function(x){x[1]}))
matxdat$spcol<-unlist(lapply(strsplit(matxdat$id, '\\.'), function(x){x[2]}))

d1<-with(matxdat[matxdat$Resample!=matxdat$spcol,c(9,10,4)], 
     structure(auc, Size = length(unique(matxdat$Resample)),
      Labels = unique(matxdat$Resample),
      Diag = F, Upper = FALSE,method = "user", class = "dist"))

d1<-(2-d1)

d2<-with(matxdat[matxdat$Resample!=matxdat$spcol,c(9,10,8)], 
         structure(TSS, Size = length(unique(matxdat$Resample)),
                   Labels = unique(matxdat$Resample),
                   Diag = F, Upper = FALSE,method = "user", class = "dist"))
d2<-(2-d2)

hc_auc<-hclust(d1, method='average')
hc_tss<-hclust(d2, method='average')

png(paste0('C:/seabirds/data/modelling/plots/',k, '_hclust.png'),width = 6, height =12 , units ="in", res =300)
par(mfrow=c(2,1))
plot(hc_auc, main='AUC-derived clustering')
plot(hc_tss, main='TSS-derived clustering')
dev.off()
par(mfrow=c(1,1))

# Visualising predicitve ability between colonies
names(all_col_tune)[names(all_col_tune)=='Resample.x']<-'spcol'
all_col_tune$sp<-NULL
all_col_tune$Resample.y<-NULL
aucz<-rbind(indiv_col_tune, data.frame(Resample='MultiCol', all_col_tune))

temp1<-aucz%>%group_by(Resample)%>%
  filter(as.character(spcol)!=as.character(Resample))%>%summarise_if(is.numeric ,sum)
aucz<-rbind(aucz,data.frame(temp1[,1], spcol='SUM', temp1[,2:8]))

aucz$Resample_auc<-factor(aucz$Resample, levels=c("MultiCol",rev(hc_auc$labels[hc_auc$order])))
aucz$spcol_auc<-factor(aucz$spcol,levels=c("SUM",rev(hc_auc$labels[hc_auc$order])))
aucz$Resample_tss<-factor(aucz$Resample, levels=c("MultiCol",rev(hc_tss$labels[hc_tss$order])))
aucz$spcol_tss<-factor(aucz$spcol,levels=c("SUM",rev(hc_tss$labels[hc_tss$order])))


aucz$auctemp<-aucz$auc
aucz[aucz$spcol=='SUM',]$auctemp<-NA
aucz$auc_bin<-cut(aucz$auctemp, c(0, 0.6, 0.7, 0.8, 0.9, 1),include.lowest =T)
aucz$auc_bin<-factor(aucz$auc_bin, levels = levels(addNA(aucz$auc_bin)),
                     labels = c("#cccccc", "#ffffb2", "#fecc5c", '#fd8d3c', '#f03b20', '#c6dbef'), exclude = NULL)

aucz$tsstemp<-aucz$TSS
aucz[aucz$spcol=='SUM',]$tsstemp<-NA
aucz$tss_bin<-cut(aucz$tsstemp, c(0, 0.2, 0.6, 1), include.lowest =T)
aucz$tss_bin<-factor(aucz$tss_bin, levels = levels(addNA(aucz$tss_bin)),
                     labels = c("#cccccc", "#fecc5c", '#f03b20', '#c6dbef'), exclude = NULL)

p_auc<-ggplot(aucz, aes(x = spcol_auc, y = Resample_auc)) + 
  geom_raster(aes(fill=auc_bin)) +scale_fill_identity()+ 
  geom_text(aes(label=round(auc, 2)), size=2)+
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+ylab('Predictions from')+xlab('Predicting to')

p_tss<-ggplot(aucz, aes(x = spcol_tss, y = Resample_tss)) + 
  geom_raster(aes(fill=tss_bin)) +scale_fill_identity()+ 
  geom_text(aes(label=round(TSS, 2)), size=2)+
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+ylab('Predictions from')+xlab('Predicting to')

# save outputs
matx_out<-rbind(matx_out, data.frame(sp=k, matxdat)) # capture hclust for plotting
aucz_out<-rbind(aucz_out, data.frame(sp=k, aucz)) # capture validation for plotting

}# close loop

# write out saved auc and hclust data
write.csv(aucz_out, 'C:/seabirds/data/mod_validation_vals.csv', quote=F, row.names=F)
write.csv(matx_out, 'C:/seabirds/data/mod_clustering_vals.csv', quote=F, row.names=F)

####~~ Individual colony model overlap ~~####
pred_list<-list.files('C:/seabirds/data/modelling/GBR_preds', full.names=T)
pred_list<-pred_list[-grep('MultiCol', pred_list)]
pred_list<-pred_list[-grep('indiv', pred_list)]
sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG', 'FRBD', 'TRBD', 'NODD', 'TERN')
for(i in sp_groups){
sp_stack<-stack(pred_list[grep(i, pred_list)])
#sp_norm<-normalized(sp_stack)
sp_sum<-calc(sp_stack, sum)

# Read in tuning results
indiv_col_tune<-read.csv(paste0('C:/seabirds/data/modelling/rf_tuning/',i, '_indiv_col_tune.csv'))
# lookup optimal vals
indiv_col_tune<-na.omit(left_join(indiv_col_tune, filter(my_hyp, sp==i),
                                  by=c('Resample', 'mtry', 'min.node.size')))

ict_med<-indiv_col_tune%>%group_by(Resample)%>%summarise(thresh=median(thresh, ra.rm=T))

for(k in 1:nlayers(sp_stack))
{
 r1<-subset(sp_stack, k)
 med_lookup<-ict_med[ict_med$Resample==gsub('_', ' ',substr(names(r1), 6,nchar(names(r1)))),]$thresh
 r1<-reclassify(r1, c(-Inf, med_lookup, 0, med_lookup,Inf,1))
 if(k==1){thresh_stack<-r1}else{thresh_stack<-stack(thresh_stack, r1)}
}

thresh_sum<-calc(thresh_stack, sum)

#sp_norm_sum<-calc(sp_norm, sum)
writeRaster(sp_sum, paste0('C:/seabirds/data/modelling/GBR_preds/', i, '_indivSUM.tif'),overwrite=T )
writeRaster(thresh_sum, paste0('C:/seabirds/data/modelling/GBR_preds/', i, '_indivSUM_class.tif'),overwrite=T)
print(i)}
 

#temp checking
dat$spcol<-paste('BRBO', dat$spcol, sep='_')
dat<-dat%>%group_by(spcol)%>%mutate(index=1:n())%>%as.data.frame()
#dat$sum_mod<-rowSums(dat[,16:31])
#dat$mean_mod<-rowMeans(dat[,16:31])
#dat<-dat%>%group_by(spcol)%>%mutate(sum_mod_norm=normalized(sum_mod))%>%as.data.frame()

d1<-dat[,c(1,2,16:32)]%>%tidyr::gather('predcol', 'pred', -spcol,-forbin,-index)

d1<-d1%>%group_by(spcol, predcol)%>%mutate(pred_norm=normalized(pred))%>%as.data.frame()

d_sum<-d1%>%filter(spcol!=predcol)%>%group_by(forbin, spcol, index)%>%
  summarise(predcol='SUM',pred=sum(pred))%>%as.data.frame()

d_sum2<-d1%>%filter(spcol!=predcol)%>%group_by(forbin, spcol, index)%>%
  summarise(predcol='SUM2',pred=sum(pred_norm))%>%as.data.frame()

d3<-rbind(d1%>%select(-pred_norm), d_sum, d_sum2)

indiv_aucz<-d3%>%group_by(spcol,predcol)%>%
  summarise(auc=as.double(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'), direction="<")$auc),
            thresh=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
indiv_aucz$TSS=indiv_aucz$sens+indiv_aucz$spec-1

dt1<-indiv_aucz%>%filter(spcol!=predcol )%>%group_by(predcol)%>%summarise(mean(auc))
dt1$auc_norm<-normalized(dt1$`mean(auc)`)
d4<-left_join(d1, dt1, by="predcol")

d_sum3<-d4%>%filter(spcol!=predcol)%>%group_by(forbin, spcol, index)%>%
  summarise(predcol='SUM3',pred=sum(pred_norm*auc_norm))%>%as.data.frame()

indiv_auczSUM<-d_sum3%>%group_by(spcol, predcol)%>%
  summarise(auc=as.double(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'), direction="<")$auc),
            thresh=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
indiv_aucz$TSS=indiv_aucz$sens+indiv_aucz$spec-1

indiv_aucz2<-rbind(indiv_aucz, indiv_auczSUM)

indiv_aucz2%>%filter(spcol!=predcol )%>%group_by(predcol)%>%summarise(mean(auc))

ggplot(indiv_aucz2, aes(x = spcol, y = predcol)) + 
       geom_raster(aes(fill=auc))+ 
       geom_text(aes(label=round(auc, 2)), size=2)+
       theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))+ylab('Predictions from')+xlab('Predicting to')
        