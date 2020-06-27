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

set.seed(24) # so that random samples and random forest models can be reproduced
normalized<-function(x){(x-min(x))/(max(x)-min(x))} # normalise (0-1) function

# read in data
sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG',
                      'FRBD', 'TRBD', 'NODD', 'TERN')

for(k in sp_groups)

dat<-read.csv(paste0('C:/seabirds/data/modelling/kernhull_pts/', k, '_kernhull.csv'))
names(dat)[names(dat)=='layer']<-'forbin'
dat$forbin<-as.factor(dat$forbin)
dat[dat$forbin==1,]$weight=0 # change from NA for no.omit

# RM sp name unless k has conflicts
if(k!='FRBD'){dat$spcol<-substr(dat$spcol, 6, nchar(as.character(dat$spcol)))}

# clean up data
# set bathy vals >0 to NA then make positive
dat[dat$bth>0,]$bth<-NA
dat$bth<-sqrt(dat$bth^2)
dat<-na.omit(dat)

dat%>%group_by(spcol)%>%summarise(nP=length(which(forbin==1)),nA=length(which(forbin==0)))

dat<-dat%>%group_by(spcol)%>%mutate(npres=length(which(forbin==1)))

# 3:1 Pa to A, sampled wothout replacement, over normalised inv coldist surface
dat<-rbind(dat%>%filter(forbin==1), dat%>%filter(forbin==0)%>%group_by(spcol)%>%
          mutate(w2=normalized(weight))%>% 
          sample_n(size=(unique(npres)*3), replace=F, weight=w2))%>%as.data.frame() 

#check ggplot(data=filter(dat, forbin==0), aes(x=Longitude, y=Latitude, colour=w2))+geom_point(size=0.4)+facet_wrap(~spcol, scales='free')
dat$npres<-NULL
dat$w2<-NULL

# niche overlap
enviro_std<-decostand(dat[,c(4:13)], method="standardize")
enviro_rda<-rda(enviro_std, scale=T)
screeplot(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1))
enviro.sites.scores<-data.frame(enviro.sites.scores,dat[,c(1,2)])
enviro.species.scores<-as.data.frame(scores(enviro_rda, display='species'))
enviro.species.scores$Predictors<-colnames(enviro_std)
enviro.species.scores$spcol='Danger'
enviro.species.scores$PC1<-enviro.species.scores$PC1/30
enviro.species.scores$PC2<-enviro.species.scores$PC2/30
ggplot(data=enviro.sites.scores, aes(x=PC1, y=PC2))+
  geom_segment(data=enviro.species.scores, aes(y=0, x=0, yend=PC2, xend=PC1), arrow=arrow(length=unit(0.3,'lines')), colour='black')+
  geom_text(data=enviro.species.scores, aes(y=PC2, x=PC1, label=Predictors), colour='black')+
  geom_density_2d(aes(colour=factor(forbin)))+facet_wrap(~spcol)+theme_bw()+
  scale_colour_manual('Area',values=c("seagreen3","coral2"), labels=c('Accessible \n habitat', 'Core \n foraging'))

# Remember nearshore front values which ==0 should be NA

# scale: normalise (0-1) each covariate per dataID
# Remember to check for outliers prior to normaliziation
#dat_norm<-dat %>% group_by(spcol) %>% mutate_if(is.numeric, normalized)%>%as.data.frame()

dat_norm<-dat

sp_store<-NULL
var_imp<-NULL
for( i in unique(dat_norm$spcol))
{
  #pr1<-table(dat_norm[dat_norm$spcol==i,]$forbin)[1]/table(dat_norm[dat_norm$ID==i,]$forbin)[2]
  cw<-ifelse(dat_norm[dat_norm$spcol==i,]$forbin==0, 1, 3)  
  
  rf1<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr+bth+slp,
              data=dat_norm[dat_norm$spcol==i,], num.trees=500, mtry=3, importance='impurity',
              probability =T)
  print(rf1)
  
  # predict to other colonies
  dat_norm$p1<-predict(rf1, data=dat_norm)$predictions[,2] # prob of foraging 0-1
  names(dat_norm)[which(names(dat_norm)=='p1')]<-i
  # predict to GBR
  #gbr_norm$p1<-predict(rf1, newdata=gbr_norm, type='prob')[,2] # prob of foraging 0-1
  #names(gbr_norm)[which(names(gbr_norm)=='p1')]<-i
  
  #varib importance
  
  var_imp<-rbind(var_imp, data.frame(spcol=i, t(ranger::importance(rf1)/max(ranger::importance(rf1)))))
  
  # rm(rf1) # save space
  
   #SPAC assessment
  wtdist<-max(dat_norm[dat_norm$spcol==i,]$weight)-100 # 100 km from col
  residz<-as.integer(as.character(dat_norm[dat_norm$spcol==i & dat_norm$weight>wtdist,]$forbin))-
    dat_norm[dat_norm$spcol==i & dat_norm$weight>wtdist,i]
  sp1<-spline.correlog(dat_norm[dat_norm$spcol==i& dat_norm$weight>wtdist,]$Longitude,
                       dat_norm[dat_norm$spcol==i& dat_norm$weight>wtdist,]$Latitude, residz,
                       na.rm=T,latlon=T,resamp=10)
  plot(sp1)
  sp_store<-rbind(sp_store, data.frame(spcol=i, Dist=sp1$boot$boot.summary$predicted$x[1,],
                       SPAC_025=sp1$boot$boot.summary$predicted$y[3,],
                       SPAC_Ave=sp1$boot$boot.summary$predicted$y[6,],
                       SPAC_95=sp1$boot$boot.summary$predicted$y[9,]))
  print(i)
  print(Sys.time()) 
}

qplot(data=sp_store, x=Dist, y=SPAC_Ave, colour=spcol, geom="line")+
  geom_ribbon(aes(ymin=SPAC_025, ymax=SPAC_95, fill=spcol),alpha=0.25)+theme_classic()

# check cols selected and size of matrix before running
ncolz<-length(unique(dat_norm$spcol))
br1<-as.data.frame(dat_norm[,c(1, 2, 16:(16+(ncolz-1)))])%>%tidyr::gather(variable, value, -spcol, -forbin)
aucz<-br1%>%group_by(spcol, variable)%>%
  summarise(auc=pROC::roc(forbin, value,direction="<")$auc,
            thresh=coords(pROC::roc(forbin, value,direction="<"),'best', best.method='closest.topleft', transpose=F)$threshold,
            sens=coords(pROC::roc(forbin, value,direction="<"),'best', best.method='closest.topleft', transpose=F)$sensitivity,
            spec=coords(pROC::roc(forbin, value,direction="<"),'best', best.method='closest.topleft', transpose=F)$specificity)
#closest.topleft works, youden shows error
aucz$TSS=aucz$sens+aucz$spec-1
# handy
qplot(data=filter(br1, spcol=='Prickly Pear'& variable=='Danger'), x=value, fill=factor(forbin), geom='histogram')


m1<-matrix(ncol=ncolz, nrow=ncolz, data =aucz$auc, dimnames=list(unique(aucz$spcol), unique(aucz$spcol)))
d1<-as.dist(1-m1)

plot(hclust(d1, method='average'))

p_auc<-ggplot(aucz, aes(x = spcol, y = variable)) + 
  geom_raster(aes(fill=auc)) + 
  geom_text(aes(label=round(auc, 2)))+
  scale_fill_gradient(low="grey",  high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+ylab('Model predictions')+xlab('Predicting to')

p_tss<-ggplot(aucz, aes(x = spcol, y = variable)) + 
  geom_raster(aes(fill=TSS)) + 
  geom_text(aes(label=round(TSS, 2)))+
  scale_fill_gradient(low="grey",  high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+ylab('Model predictions')+xlab('Predicting to')


br2<-as.data.frame(dat_norm[,c(1, 2, 14,15,16:(16+(ncolz-1)))])%>%
  tidyr::gather(variable, value, -spcol, -forbin, -Longitude, -Latitude)
br2<-left_join(br2, aucz[,c(1,2,4)], by=c('spcol', 'variable'))

ggplot(data=br2, aes(x=Longitude, y=Latitude))+
  geom_point(aes(colour=value))+geom_point(data=br2[br2$value>unique(br2$thresh),], shape=1, colour=2)+
  facet_grid(spcol~variable)

#caret approach

#test_dat<-as.data.frame(dat_norm[dat_norm$ID %in% unique(dat_norm$ID)[c(2,4,5)],])

#try subsample data to something that will actually run
dat_norm%>%group_by(ID)%>%summarise(n())
test_dat<-dat_norm%>%group_by(ID)%>%
  slice(1:4500)%>%as.data.frame()  # grab first 4500 rows

# test_dat has to be data.frame otherwise caret 'y' classerror

# make forbin a factor with labels attached rather than 0, 1
test_dat$forbin<-as.factor(test_dat$forbin)
levels(test_dat$forbin)<-c("transit", "forage")

# do somne variable vis
td <- test_dat %>% dplyr::select(-embc, -trip_id)%>%gather(variable, value, -ID, -forbin)
ggplot(data=td, aes(x=value, colour=forbin))+
  geom_density()+facet_grid(ID~variable, scales='free')


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

# PCA
enviro_std<-decostand(test_dat[,c(3:6, 8, 9, 11, 12)], method="standardize")
enviro_rda<-rda(enviro_std, scale=T)
screeplot(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1))
enviro.sites.scores<-data.frame(enviro.sites.scores,test_dat[,c(1,13)])
ggplot(data=enviro.sites.scores, aes(x=PC1, y=PC2))+
  geom_point(aes(colour=forbin), shape=1, alpha=0.5)+facet_wrap(~ID)
