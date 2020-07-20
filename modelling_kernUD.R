# Modelling (kernelUD approach)

library(dplyr)
library(tidyr)
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

#read in GBR ocean data
gbr<-read.csv('C:/seabirds/data/pred_area_modelready_2km.csv')
#rename
names(gbr)<-c('x', 'y', 'sst', 'sst_sd', 'chl', 'chl_sd',
              'mfr', 'mfr_sd', 'pfr', 'pfr_sd', 'bth', 'slp')
# and clean
gbr[gbr$bth>0,]$bth<-0
gbr$bth<-sqrt(gbr$bth^2)

# read in 2km rasterize template
templ<-raster('C:/seabirds/data/GIS/pred_area_ras_template2km.tif')

# read in data
sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

for(k in sp_groups)
{
dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts_sample/', k, '_kernhull_sample.csv'))
dat$X<-NULL

# niche overlap
enviro_std<-decostand(dat[,c(4:13)], method="standardize")
enviro_rda<-rda(enviro_std, scale=T)
screeplot(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1))
enviro.sites.scores<-data.frame(enviro.sites.scores,dat[,c(1,2)])
enviro.species.scores<-as.data.frame(scores(enviro_rda, display='species'))
enviro.species.scores$Predictors<-colnames(enviro_std)
enviro.species.scores$spcol='X-Predictors'
enviro.species.scores$PC1<-enviro.species.scores$PC1/30
enviro.species.scores$PC2<-enviro.species.scores$PC2/30

pniche<-ggplot(data=enviro.sites.scores, aes(x=PC1, y=PC2))+
  geom_segment(data=enviro.species.scores, aes(y=0, x=0, yend=PC2, xend=PC1), arrow=arrow(length=unit(0.3,'lines')), colour='black')+
  geom_text(data=enviro.species.scores, aes(y=PC2, x=PC1, label=Predictors), colour='black', size=2)+
  geom_density_2d(aes(colour=forbin), alpha=0.6)+facet_wrap(~spcol)+theme_bw()+
  scale_colour_manual('Area',values=c("coral2", "seagreen3"), labels=c('Core \n foraging', 'Accessible \n habitat'))

# Do hull and order sites by dendrogram or geometric distance

#png(paste0('C:/seabirds/data/modelling/plots/',k, '_niche.png'),width = 6, height =6 , units ="in", res =600)
#print(pniche)
#dev.off()

# Remember nearshore front values which ==0 should be NA

# scale: normalise (0-1) each covariate per dataID
# Remember to check for outliers prior to normaliziation
#dat<-dat %>% group_by(spcol) %>% mutate_if(is.numeric, normalized)%>%as.data.frame()

sp_store<-NULL
var_imp<-NULL
for( i in unique(dat$spcol))
{
  #pr1<-table(dat[dat$spcol==i,]$forbin)[1]/table(dat[dat$ID==i,]$forbin)[2]

  rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp,
              data=dat[dat$spcol==i,], num.trees=500, 
              mtry=my.mtry,  splitrule = "gini", min.node.size =my.mns,  importance='impurity',
              probability =T)
  print(rf1)
  
  # predict to other colonies
  dat$p1<-predict(rf1, data=dat)$predictions[,2] # prob of foraging 0-1
  names(dat)[which(names(dat)=='p1')]<-i
  # predict to GBR
  #gbr_norm$p1<-predict(rf1, newdata=gbr_norm, type='prob')[,2] # prob of foraging 0-1
  #names(gbr_norm)[which(names(gbr_norm)=='p1')]<-i
  
  #varib importance
  
  var_imp<-rbind(var_imp, data.frame(spcol=i, t(ranger::importance(rf1)/max(ranger::importance(rf1)))))
  
  # predict to GBR
  gbr$p1<-predict(rf1, data=gbr)$predictions[,2] # prob of foraging 0-1
  names(gbr)[which(names(gbr)=='p1')]<-paste(k, i, sep='_')
  
  # rm(rf1) # save space
  
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
}

# make full model and predict to GBR

rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp ,
            data=dat, num.trees=500, mtry=3, importance='impurity',
            probability =T)
gbr$MultiCol<-predict(rf1, data=gbr)$predictions[,2]
names(gbr)[which(names(gbr)=='MultiCol')]<-paste(k, 'MultiCol', sep='_')

spdf<-SpatialPointsDataFrame(SpatialPoints(gbr[,1:2], proj4string = CRS(projection(templ))),
                             data=gbr[,grep(k, names(gbr))])

for(j in 1:length(spdf@data))
{
  p1<-rasterize(spdf, templ, field=names(spdf@data)[j])
  writeRaster(p1, paste0('C:/seabirds/data/modelling/GBR_preds/', gsub(' ', '_', names(spdf@data)[j]), '.tif'), overwrite=T)
  #if(i==1){ps<-p1}else{ps<-stack(ps, p1)}
}

#print(k)
#}


#qplot(data=sp_store, x=Dist, y=SPAC_Ave, colour=spcol, geom="line")+
#  geom_ribbon(aes(ymin=SPAC_025, ymax=SPAC_95, fill=spcol),alpha=0.25)+theme_classic()

#write.csv(sp_store, paste0('C:/seabirds/data/modelling/SPAC/',k,'_spac.csv'), quote=F, row.names = F)
write.csv(var_imp, paste0('C:/seabirds/data/modelling/var_imp/',k,'_var_imp.csv'), quote=F, row.names = F)

# check cols selected and size of matrix before running
ncolz<-length(unique(dat$spcol))
br1<-as.data.frame(dat[,c(1, 2, 16:(16+(ncolz-1)))])%>%tidyr::gather(variable, value, -spcol, -forbin)
aucz<-br1%>%group_by(spcol, variable)%>%
  summarise(auc=as.double(pROC::roc(forbin, value,direction="<")$auc),
            thresh=coords(pROC::roc(forbin, value,direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(forbin, value,direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(forbin, value,direction="<"),'best', best.method='youden', transpose=F)$specificity[1])

#closest.topleft works or youden 
# calc TSS
aucz$TSS=aucz$sens+aucz$spec-1
# handy
#qplot(data=filter(br1, spcol=='Prickly Pear'& variable=='Danger'), x=value, fill=forbin, geom='histogram')

m1<-matrix(ncol=ncolz, nrow=ncolz, data =aucz$auc, dimnames=list(unique(aucz$spcol), unique(aucz$spcol)))
d1<-as.dist(1-m1)

m1<-matrix(ncol=ncolz, nrow=ncolz, data =aucz$TSS, dimnames=list(unique(aucz$spcol), unique(aucz$spcol)))
d2<-as.dist(1-m1)

png(paste0('C:/seabirds/data/modelling/plots/',k, '_hclust.png'),width = 6, height =12 , units ="in", res =300)
par(mfrow=c(2,1))
plot(hclust(d1, method='average'), main='AUC-derived clustering')
plot(hclust(d2, method='average'), main='TSS-derived clustering')
dev.off()
par(mfrow=c(1,1))


# caret LGOCV combi-spcol model
folds <- groupKFold(dat$spcol)

train_control <- trainControl( method="LGOCV",index=folds,
                               classProbs = TRUE, savePredictions = TRUE,
                               summaryFunction = twoClassSummary)

tunegrid <- expand.grid(mtry=3,  splitrule = "gini", min.node.size = 10)

rf_default <- caret::train(x=dat[,c('sst','sst_sd','chl','chl_sd','mfr_sd',
                                         'pfr','bth','slp')],
                           y=dat[,'forbin'], method="ranger", num.trees=500, metric='ROC', 
                           tuneGrid=tunegrid, trControl=train_control)
print(rf_default)

rf_default$resample

head(rf_default$pred)

caret_aucz<-rf_default$pred%>%group_by(Resample)%>%
  summarise(auc=as.double(pROC::roc(obs, Core, direction="<")$auc),
            thresh=coords(pROC::roc(obs, Core,direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(obs, Core,direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(obs, Core,direction="<"),'best', best.method='youden', transpose=F)$specificity[1])

caret_aucz<-data.frame(spcol=paste(unlist(lapply(folds, function(x){unique(dat$spcol) [-which(unique(dat$spcol) %in% unique(dat[x,]$spcol))][[1]]}))),
                       variable='MultiCol', caret_aucz[,2:5])
caret_aucz$TSS=caret_aucz$sens+caret_aucz$spec-1

aucz<-rbind(as.data.frame(aucz), as.data.frame(caret_aucz))
aucz<-rbind(aucz, data.frame(spcol='SUM', aucz%>%group_by(variable)%>%
                               filter(as.character(spcol)!=as.character(variable))%>%summarise_if(is.numeric ,sum)))
aucz$variable<-factor(aucz$variable)
aucz$variable<-relevel(aucz$variable, ref = "MultiCol")
aucz$spcol<-factor(aucz$spcol)
aucz$spcol<-relevel(aucz$spcol, ref = "SUM")

p_auc<-ggplot(aucz, aes(x = spcol, y = variable)) + 
  geom_raster(aes(fill=auc)) + 
  geom_text(aes(label=round(auc, 2)), size=2)+
  scale_fill_gradient(limits=c(0.5, 1), low="grey",  high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+ylab('Model predictions')+xlab('Predicting to')

p_tss<-ggplot(aucz, aes(x = spcol, y = variable)) + 
  geom_raster(aes(fill=TSS)) + 
  geom_text(aes(label=round(TSS, 2)), size=2)+
  scale_fill_gradient(limits=c(0,1),low="grey",  high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+ylab('Model predictions')+xlab('Predicting to')

png(paste0('C:/seabirds/data/modelling/plots/',k, '_pred_acc.png'),width = 6, height =12 , units ="in", res =600)
grid.arrange(p_auc, p_tss)
dev.off()

}# close loop

br2<-as.data.frame(dat[,c(1, 2, 14,15,16:(16+(ncolz-1)))])%>%
  tidyr::gather(variable, value, -spcol, -forbin, -Longitude, -Latitude)
br2<-left_join(br2, aucz[,c(1,2,4)], by=c('spcol', 'variable'))
br2$pred_core<-ifelse(br2$value>=br2$thresh, 1, 0)

ggplot(data=br2[br2$spcol=='Swains',], aes(x=Longitude, y=Latitude))+
  geom_point(aes(colour=value), size=0.5)+
  geom_point(data=br2[br2$spcol=='Swains' & br2$pred_core==1,], shape=1, colour='green', size=0.5)+
  scale_colour_gradientn(colours = heat.colors(10))+
  facet_wrap(~variable, scales='free')+theme_bw()

#caret approach

#test_dat<-as.data.frame(dat[dat$ID %in% unique(dat$ID)[c(2,4,5)],])

#try subsample data to something that will actually run
dat%>%group_by(ID)%>%summarise(n())
test_dat<-dat%>%group_by(ID)%>%
  slice(1:4500)%>%as.data.frame()  # grab first 4500 rows

# test_dat has to be data.frame otherwise caret 'y' classerror

# make forbin a factor with labels attached rather than 0, 1
test_dat$forbin<-as.factor(test_dat$forbin)
levels(test_dat$forbin)<-c("transit", "forage")



test_dat$ID<-factor(test_dat$ID)
folds <- groupKFold(test_dat$ID)

train_control <- trainControl( method="LGOCV",
                               index=folds,
                               classProbs = TRUE,
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary)

tunegrid <- expand.grid(mtry=3,  
                        splitrule = "gini",
                        min.node.size = 10)

### DONT USE FORMULA INTERFACE!

rf_default <- caret::train(x=test_dat[,c('sst_mn','sst_sd','chl_mn','chl_sd','mfront_sd',
                                         'pfront_mn','ex_bathy','ex_slope')],
                           y=test_dat[,12], method="ranger", num.trees=500, metric='ROC', 
                           tuneGrid=tunegrid, trControl=train_control)
print(rf_default)

rf_default$resample

head(rf_default$pred)

rf_default$pred%>%group_by(Resample)%>%summarise(auc=pROC::roc(obs, forage)$auc,
                                                 sens=pROC::roc(obs, forage)$sensitivities[which(pROC::roc(obs, forage)$thresholds>0.5 & pROC::roc(obs, forage)$thresholds<0.51)[1]],
                                                 spec=pROC::roc(obs, forage)$specificities[which(pROC::roc(obs, forage)$thresholds>0.5 & pROC::roc(obs, forage)$thresholds<0.51)[1]])

# Check why thereis huge difference between sensitivity and specificity, 
# basically the model isn't predicting presences
table(rf_default$pred$pred); table(rf_default$pred$obs)


# OLD
# do somne variable vis
td <- test_dat %>% dplyr::select(-embc, -trip_id)%>%gather(variable, value, -ID, -forbin)
ggplot(data=td, aes(x=value, colour=forbin))+
  geom_density()+facet_grid(ID~variable, scales='free')

