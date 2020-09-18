# Modelling (kernelUD approach)

library(dplyr)
library(ggplot2)
library(ranger)
library(pROC)
library(raster)

set.seed(24) # so that random samples and random forest models can be reproduced
normalized<-function(x){(x-min(x))/(max(x)-min(x))} # normalise (0-1) function

####~~ Random sampling of psuedo absences fixed for reproducability amd data compression ~~~####
#sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG','FRBD', 'TRBD', 'NODD', 'TERN')
#for(k in sp_groups)
#{
#  for(m in 2:5){
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
#  dat%>%group_by(spcol)%>%summarise(nP=length(which(forbin=='Core')),nA=length(which(forbin=="PsuedoA"))) # don't need
#    dat<-dat%>%group_by(spcol)%>%mutate(npres=length(which(forbin=='Core')))

    # 3:1 Pa to A, sampled wothout replacement, over normalised inv coldist surface
#  dat<-rbind(dat%>%filter(forbin=='Core'), dat%>%filter(forbin=="PsuedoA")%>%group_by(spcol)%>%
#               mutate(w2=normalized(weight))%>% 
#               sample_n(size=(unique(npres)*3), replace=F, weight=w2))%>%as.data.frame() 
#  dat$npres<-NULL
#  dat$w2<-NULL
#  write.csv(dat,paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', k, '_kernhull_sample', m, '.csv'), row.names=F, quote=F)
#  print(k)
#  }
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
gbr<-read.csv('C:/seabirds/data/pred_area_large_modelready_2km.csv')
names(gbr)<-c('x', 'y', 'sst', 'sst_sd', 'chl', 'chl_sd',
              'mfr', 'mfr_sd', 'pfr', 'pfr_sd', 'bth', 'slp')
gbr[gbr$bth>0,]$bth<-0 
gbr$bth<-sqrt(gbr$bth^2)# Remember nearshore front values which ==0 should be NA

# read in 2km rasterize template
templ<-raster('C:/seabirds/data/GIS/pred_area_large_ras_template.tif')

# read in data
sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

matx_out<-NULL # capture hclust for plotting
aucz_out<-NULL # capture validation for plotting

for(k in sp_groups)
{
  # Using the first resampled dataset as the baseline for AUC validation 
  dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', k, '_kernhull_sample1.csv'))
  dat$X<-NULL

  # lookup optimal vals
  my.hyp.sp<-filter(my_hyp, sp==k)
  my.hyp.sp$Resample<-as.character(my.hyp.sp$Resample)
  
  if(k=='RFBO'){dat<-filter(dat, spcol!='Christmas')}
  if(k=='SOTE'){dat$spcol<-as.character(dat$spcol)
  dat[dat$spcol=='chick',]$spcol<-'Rat'}
  if(k=='WTST'| k=='WTLG'){
    gtemp<-read.csv('C:/seabirds/data/pred_area_large_modelready_2km_WTSHsummer.csv')
    gbr2<-gbr
    gbr2$sst<-gtemp$ex_sst
    gbr2$chl<-gtemp$ex_chl
    gbr2$mfr<-gtemp$ex_mfr
    gbr2$pfr<-gtemp$ex_pfr
    rm(gtemp)}else{gbr2<-gbr}

dat$MultiCol<-0 # filled for each LGOCV in loop  
rep_preds<-data.frame(ID=1:nrow(dat))  
rep_GBRpreds<-data.frame(ID=1:nrow(gbr2))
var_imp<-NULL
for(h in 1:5) # loop through 5 resamples of PsuedoAbs data and evrage mode predictions
{
  if(h==1){repdat<-dat}else{
    repdat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', k, '_kernhull_sample', h,'.csv'))
  } 
  for( i in unique(repdat$spcol))
  {
      
    #predict colony model to other colonies
    rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp,
                  data=repdat[repdat$spcol==i,], num.trees=500, 
                  mtry=filter(my.hyp.sp, Resample=='MultiCol')$mtry, 
                  min.node.size=filter(my.hyp.sp, Resample=='MultiCol')$min.node.size,
                  splitrule = "gini",  importance='impurity',probability =T, seed=24)  
      
    # predict to other colonies 
    rep_preds$p1<-predict(rf1, data=dat)$predictions[,1] # prob of foraging 0-1
    names(rep_preds)[which(names(rep_preds)=='p1')]<-paste(i, h, sep='_')
    
    # predict to GBR
    rep_GBRpreds$p1<-predict(rf1, data=gbr2)$predictions[,1] # prob of foraging 0-1
    names(rep_GBRpreds)[which(names(rep_GBRpreds)=='p1')]<-paste(i, h, sep='_')
   
    #varib importance
    var_imp<-rbind(var_imp, data.frame(spcol=i, rep=h, t(ranger::importance(rf1)/max(ranger::importance(rf1)))))
    
    # make leave-one-out multicolony model
    rf2<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp ,
                data=repdat[repdat$spcol!=i,], num.trees=500,  
                mtry=filter(my.hyp.sp, Resample=='MultiCol')$mtry, 
                min.node.size=filter(my.hyp.sp, Resample=='MultiCol')$min.node.size,
                splitrule = "gini",  importance='impurity',probability =T, seed=24)
    
    # predict to other colonies 
    rep_preds$MultiCol<-predict(rf2, data=dat)$predictions[,1] # prob of foraging 0-1
    names(rep_preds)[which(names(rep_preds)=='MultiCol')]<-paste('MultiCol',i, h, sep='_')
    
    #averageing 1:5 predictions for model validation and GBR pred
    if(h==5)
        {
        # mod val
        dat$p3<-rowMeans(rep_preds[paste(i, 1:5, sep='_')])
        names(dat)[which(names(dat)=='p3')]<-i
        dat[dat$spcol==i,]$MultiCol<-rowMeans(rep_preds[paste('MultiCol',i,1:5,sep='_')])[which(dat$spcol==i)]
        #spcol gbr pred
        gbr2$p3<-rowMeans(rep_GBRpreds[paste(i, 1:5, sep='_')])
        names(gbr2)[which(names(gbr2)=='p3')]<-i
        }
    print(i)
    rm(rf1) # save space
    rm(rf2)
    print(head(rep_preds))  
  }

# make full (all cols) multicolony model for prediction to Reef
  rf3<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp ,
              data=repdat, num.trees=500,  
              mtry=filter(my.hyp.sp, Resample=='MultiCol')$mtry, 
              min.node.size=filter(my.hyp.sp, Resample=='MultiCol')$min.node.size,
              splitrule = "gini",  importance='impurity',probability =T, seed=24)
 
  #predict to GBR
  rep_GBRpreds$MultiCol<-predict(rf3, data=gbr2)$predictions[,1]
  names(rep_GBRpreds)[which(names(rep_GBRpreds)=='MultiCol')]<-paste('MultiCol', h, sep='_')
  
  if(h==5) # ave 1:5 multicol model preds to GBR
  {gbr2$MultiCol<-rowMeans(rep_GBRpreds[paste('MultiCol', 1:5, sep='_')])}
    
  rm(rf2)
  print(h)
  print(Sys.time())
}

#Average 1:5 predictions
for(m in unique(dat$spcol))
indiv_aucz<-tidyr::gather(rep_preds, Resample, pred, -ID)

# Write out spatial predictions
spdf<-SpatialPointsDataFrame(SpatialPoints(gbr2[,1:2], proj4string = CRS(projection(templ))),
                             data=gbr2[,grep(k, names(gbr2))])

for(j in 1:length(spdf@data))
{
  p1<-rasterize(spdf, templ, field=names(spdf@data)[j])
  writeRaster(p1, paste0('C:/seabirds/data/modelling/GBR_preds/', gsub(' ', '_', names(spdf@data)[j]), '.tif'), overwrite=T)
  #if(i==1){ps<-p1}else{ps<-stack(ps, p1)}
}

#write.csv(sp_store, paste0('C:/seabirds/data/modelling/SPAC/',k,'_spac.csv'), quote=F, row.names = F)
write.csv(var_imp, paste0('C:/seabirds/data/modelling/var_imp/',k,'_var_imp.csv'), quote=F, row.names = F)

#calc AUC and TSS

indiv_aucz<-tidyr::gather(dat[,c(1,2,16:ncol(dat))], Resample, pred, -forbin, -spcol)%>%
  group_by(Resample, spcol)%>%
  summarise(auc=as.double(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'), direction="<")$auc),
            thresh=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
indiv_aucz$TSS=indiv_aucz$sens+indiv_aucz$spec-1

# Clustering sites based on predictive ability

indiv_aucz$id<-apply(as.data.frame(indiv_aucz), 1, FUN=function(x){paste(sort(c(as.character(x[1]), as.character(x[2]))), sep='.', collapse='.')})
matxdat<-indiv_aucz%>%filter(Resample!='MultiCol')%>%group_by(id)%>%summarise_if(is.numeric, sum)
indiv_aucz$id<-NULL
matxdat$Resample<-unlist(lapply(strsplit(matxdat$id, '\\.'), function(x){x[1]}))
matxdat$spcol<-unlist(lapply(strsplit(matxdat$id, '\\.'), function(x){x[2]}))

# calc AUC weighted ensemble of colony models

dat<-dat%>%group_by(spcol)%>%mutate(index=1:n())%>%as.data.frame()

d1<-dat[,c(1,2,16:ncol(dat))]%>%tidyr::gather('Resample', 'pred', -spcol,-forbin,-index)

d1<-d1%>%filter(Resample!='MultiCol')%>%group_by(spcol, Resample)%>%mutate(pred_norm=normalized(pred))%>%as.data.frame()

auc_mn<-indiv_aucz%>%filter(spcol!=Resample)%>%group_by(Resample)%>%
  summarise_if(is.numeric ,mean)
auc_mn$auc_norm<-(normalized(auc_mn$auc)+1)
auc_mn$auc_int<-auc_mn$auc*10
auc_mn$auc_cube<-(auc_mn$auc+1)^3
d2<-left_join(d1, auc_mn[,c(1,2,7,8, 9)], by="Resample")

ensem<-d2%>%filter(spcol!=Resample)%>%group_by(forbin, spcol, index)%>%
  summarise(Resample='EnsembleNrm',pred2=sum(pred_norm*auc_norm), pred3=sum(pred*auc_norm))%>%as.data.frame()

indiv_auczENS<-ensem%>%group_by(spcol, Resample)%>%
  summarise(auc=as.double(pROC::roc(forbin, pred2, levels=c('PsuedoA', 'Core'), direction="<")$auc),
            thresh=coords(pROC::roc(forbin, pred2, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(forbin, pred2, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(forbin, pred2, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
indiv_auczENS$TSS=indiv_auczENS$sens+indiv_auczENS$spec-1

indiv_auczENS2<-ensem%>%group_by(spcol, Resample)%>%
  summarise(auc=as.double(pROC::roc(forbin, pred3, levels=c('PsuedoA', 'Core'), direction="<")$auc),
            thresh=coords(pROC::roc(forbin, pred3, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(forbin, pred3, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(forbin, pred3, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
indiv_auczENS2$TSS=indiv_auczENS2$sens+indiv_auczENS2$spec-1
indiv_auczENS2$Resample<-'EnsembleRaw'

aucz<-rbind(indiv_aucz, indiv_auczENS, indiv_auczENS2)

# remake and bind in MEAN 
auc_mn<-aucz%>%filter(spcol!=Resample)%>%group_by(Resample)%>%
  summarise_if(is.numeric ,mean)
aucz<-rbind(as.data.frame(aucz),data.frame(auc_mn[,1], spcol='MEAN', auc_mn[,2:6]))

aucz$auctemp<-aucz$auc
aucz[aucz$spcol=='MEAN',]$auctemp<-NA
aucz$auc_bin<-cut(aucz$auctemp, c(0, 0.6, 0.7, 0.8, 0.9, 1),include.lowest =T)
aucz$auc_bin<-factor(aucz$auc_bin, levels = levels(addNA(aucz$auc_bin)),
                     labels = c("#cccccc", "#ffffb2", "#fecc5c", '#fd8d3c', '#f03b20', '#c6dbef'), exclude = NULL)

aucz$tsstemp<-aucz$TSS
aucz[aucz$spcol=='MEAN',]$tsstemp<-NA
aucz$tss_bin<-cut(aucz$tsstemp, c(0, 0.2, 0.6, 1), include.lowest =T)
aucz$tss_bin<-factor(aucz$tss_bin, levels = levels(addNA(aucz$tss_bin)),
                     labels = c("#cccccc", "#fecc5c", '#f03b20', '#c6dbef'), exclude = NULL)

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
 
####~~~*~~~####

####~~~~ Internal Block Validation for GBR-local tracking ~~~~####

# read in data
sp_groups <- c('BRBO', 'MABO', 'WTST', 'WTLG','NODD', 'TERN')

intval_out<-NULL
for(k in sp_groups)
{
  
  dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', k, '_kernhull_sample.csv'))
  dat$X<-NULL
  
  # lookup optimal vals
  my.hyp.sp<-filter(my_hyp, sp==k)
  
  if(k=='BRBO'){dcol<-c('Swains', 'Raine')}
  if(k=='MABO'){dcol<-'Swains'}
  if(k=='WTST'){dcol<-'Heron'}
  if(k=='WTLG'){dcol<-'Heron'}
  if(k=='NODD'){dcol<-'Heron'}
  
  for( i in dcol)
  {
    # Internval cross validation for GBR models
    coltemp<-dat[dat$spcol==i,]
    coltemp$clust<-kmeans(coltemp[,14:15], 4)$cluster
 
    plot(Latitude~Longitude, coltemp, col=coltemp$clust)
    points(Latitude~Longitude, coltemp[coltemp$forbin=='Core',], pch=16, cex=0.4, col=6)
    
    folds2 <- groupKFold(coltemp$clust)
    tunegrid <- expand.grid(mtry=c(2:6),  splitrule = "gini", min.node.size = c(5,10,20,50))
    
    train_control <- trainControl( method="LGOCV",index=folds2,
                                   classProbs = TRUE, savePredictions = TRUE,
                                   summaryFunction = twoClassSummary, verboseIter = TRUE)
    
    rf_intcol <- caret::train(x=coltemp[,c('sst','sst_sd','chl','chl_sd','mfr_sd', 'pfr_sd',
                                       'mfr','pfr','bth','slp')],
                              y=coltemp[,'forbin'], method="ranger", num.trees=500, metric='ROC', 
                              tuneGrid=tunegrid, trControl=train_control, verbose=T)
    
    intval_out<-rbind(intval_out,
                      data.frame(sp=k, col=i, rf_intcol$results[which.max(rf_intcol$results$ROC),]))
  }
print(intval_out)
}
write.csv(intval_out, 'C:/seabirds/data/GBR_tracking_representivity.csv', quote=F, row.names=F)
####~~~~*~~~~####
  
  
