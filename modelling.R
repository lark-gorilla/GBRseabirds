# Modelling code

library(dplyr)
library(tidyr)
library(ggplot2)
#library(GGally)
library(ranger)
library(caret)
library(pROC)
library(vegan)
library(raster)

#set.seed(7) really important to set.seed with RF??

# read in data
brbo<-read.csv('C:/seabirds/sourced_data/tracking_data/tracking_master_forage_extract.csv')
brbo<-brbo[brbo$sp=='BRBO',]

hulls<-read.csv('C:/seabirds/data/BRBO_hulls_modelready.csv')

gbr<-read.csv('C:/seabirds/data/pred_area_modelready.csv')


# clean up data

# Collapse ID to just remove owner and species, this combines all Soanes Dog
# and Swains cols - check for differences in interval
brbo$ID<-do.call(c, lapply(strsplit(as.character(brbo$ID), '_'), function(x)x[3]))

# remove sitting behaviour and error pts

brbo<-filter(brbo, embc %in% c('commuting', 'relocating', 'foraging') )
brbo$forbin<-1
brbo[brbo$embc=='commuting',]$forbin<-0

# remove NA pts (from chl, sst and front land mask mainly )
brbo$trip_id<-as.character(brbo$trip_id)
# set bathy vals >0 to NA then make positive
brbo[brbo$bth>0,]$bth<-NA
brbo$bth<-sqrt(brbo$bth^2)
brbo<-na.omit(brbo)

# clean pred area too
gbr[gbr$ex_bathy>0,]$ex_bathy<-0
gbr$ex_bathy<-sqrt(gbr$ex_bathy^2)
# quick fix to set area above 3000 to 3000
gbr[gbr$ex_bathy>3000,]$ex_bathy<-3000

# hull abs
# Collapse ID to just remove owner and species, this combines all Soanes Dog
# and Swains cols - check for differences in interval
#hulls$ID<-do.call(c, lapply(strsplit(as.character(hulls$ID), '_'), function(x)x[3]))
#hulls$ID<-substr(hulls$ID, 1, (nchar(hulls$ID)-4)) # rm .csv
#n_pres<-brbo%>%filter(forbin==1)%>%group_by(ID)%>%summarise(n())
#hull_abs<-hulls%>%group_by(ID)%>%
#  sample_n(size=filter(n_pres, ID==unique(hulls$ID))$`n()`, replace=F)

# trial to rm any duplicate points within presences and absences points (that have identical env signature)
brbo<-brbo%>%mutate(index=paste0(sst,sst_sd,chl,chl_sd,mfr_sd,pfr,bth, slp, forbin))%>%
  group_by(ID)%>%filter(!duplicated(index))%>%as.data.frame()


# trial to rm transit points that have identical env signature to foraging PER trip
brbo<-brbo%>%mutate(index=paste0(sst,sst_sd,chl,chl_sd,mfr_sd,pfr,bth, slp))%>%
  group_by(ID)%>%filter(!duplicated(index) | forbin==1)%>%as.data.frame()
# run with ID for now # loses 8k absences

# niche overlap
enviro_std<-decostand(brbo[,c(3:6, 8, 9, 11, 12)], method="standardize")
enviro_rda<-rda(enviro_std, scale=T)
screeplot(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1))
enviro.sites.scores<-data.frame(enviro.sites.scores,brbo[,c(1,13)])
ggplot(data=enviro.sites.scores, aes(x=PC1, y=PC2))+
  geom_density_2d(aes(colour=factor(forbin)))+facet_wrap(~ID)+theme_bw()


# Remember nearshore front values which ==0 should be NA
# Remember some bad tracks that should be filtered prior to modelling
# rememeber could remove bathy>0 vals

# scale: normalise (0-1) each covariate per dataID
# Remember to check for outliers prior to normaliziation
normalized<-function(x){(x-min(x))/(max(x)-min(x))}
brbo_norm<-brbo %>% group_by(ID) %>% mutate_if(is.numeric, normalized)%>%as.data.frame()
brbo_norm$forbin<-brbo$forbin
brbo$trip_id<-brbo$trip_id


for( i in unique(brbo_norm$ID))
{

#pr1<-table(brbo_norm[brbo_norm$ID==i,]$forbin)[1]/table(brbo_norm[brbo_norm$ID==i,]$forbin)[2]
#cw<-ifelse(brbo_norm[brbo_norm$ID==i,]$forbin==0, 1, pr1)  
  
rf1<-ranger(factor(forbin)~sst+sst_sd+chl+chl_sd+mfr_sd+pfr+bth+slp,
data=brbo_norm[brbo_norm$ID==i,], num.trees=500, mtry=3, importance='impurity',
probability =T)
print(rf1)
print(ranger::importance(rf1))
# predict to other colonies
brbo_norm$p1<-predict(rf1, data=brbo_norm)$predictions[,2] # prob of foraging 0-1
names(brbo_norm)[which(names(brbo_norm)=='p1')]<-i
# predict to GBR
#gbr_norm$p1<-predict(rf1, newdata=gbr_norm, type='prob')[,2] # prob of foraging 0-1
#names(gbr_norm)[which(names(gbr_norm)=='p1')]<-i
print(i)
print(Sys.time())
}

br1<-as.data.frame(brbo_norm[,c(1, 19, 21:36)])%>%gather(variable, value, -ID, -forbin)
aucz<-br1%>%group_by(ID, variable)%>%
  summarise(auc=pROC::roc(forbin, value,direction="<")$auc,
            thresh=coords(pROC::roc(forbin, value,direction="<"),0.5, transpose=F)$threshold,
            sens=coords(pROC::roc(forbin, value,direction="<"),0.5, transpose=F)$sensitivity,
            spec=coords(pROC::roc(forbin, value,direction="<"),0.5, transpose=F)$specificity)
# handy
qplot(data=filter(br1, ID=='Prickly Pear'& variable=='Danger'), x=value, fill=factor(forbin), geom='histogram')


m1<-matrix(ncol=16, nrow=16, data =aucz$auc, dimnames=list(unique(aucz$ID), unique(aucz$ID)))
d1<-as.dist(1-m1)

plot(hclust(d1, method='average'))

ggplot(aucz, aes(x = ID, y = variable)) + 
  geom_raster(aes(fill=auc)) + 
  geom_text(aes(label=round(auc, 2)))+
  scale_fill_gradient(low="grey",  high="red") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+ylab('Model predictions')+xlab('Predicting to')

# normalize GBR prediction area too
gbr_norm<-data.frame(gbr[,1:2], gbr[,3:12] %>% mutate_all(normalized))
# do bathy with custom norm to exclude set new min
#gbr_norm$ex_bathy<-(gbr$ex_bathy- -3000)/(max(x)-min(x))} not done but handy

# check correlation
#ggpairs(brbo_norm[,3:14]) # takes ages

cor(brbo_norm[,4:14])
# chl and bathy cor = 0.54 but keep
# fronts metrics cor. Keep pfront_mn and mfront_sd @ 0.47

### DONT USE FORMULA INTERFACE!

 

# Visualise GBR raster preds

templ<-raster('C:/seabirds/data/GIS/pred_area_ras_template.tif')

spdf<-SpatialPointsDataFrame(SpatialPoints(gbr_norm[,1:2], proj4string = CRS(projection(templ))),
                             data=gbr_norm[,13:length(gbr_norm)])

for(i in 1:length(spdf@data))
{
p1<-rasterize(spdf, templ, field=names(spdf@data)[i])
writeRaster(p1, paste0('C:/seabirds/data/modelling/GBR_preds/BRBO/',gsub(' ', '_', names(spdf@data)[i]), '.tif'), overwrite=T)
#if(i==1){ps<-p1}else{ps<-stack(ps, p1)}
}

#caret approach

#test_dat<-as.data.frame(brbo_norm[brbo_norm$ID %in% unique(brbo_norm$ID)[c(2,4,5)],])

#try subsample data to something that will actually run
brbo_norm%>%group_by(ID)%>%summarise(n())
test_dat<-brbo_norm%>%group_by(ID)%>%
  slice(1:4500)%>%as.data.frame()  # grab first 4500 rows


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

rf_default <- caret::train(x=test_dat[,c('ColDist','sst','sst_sd','chl', 'chl_sd',
                                         'mfr_sd', 'pfr','bth','slp')],
                           y=test_dat[,'forbin'], method="ranger", num.trees=500, metric='ROC', 
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

# full dset

# raw
car_brbo<-brbo
car_brbo$forbin<-as.factor(car_brbo$forbin)
levels(car_brbo$forbin)<-c("transit", "forage")
car_brbo$ID<-factor(car_brbo$ID)
folds <- groupKFold(car_brbo$ID)

# norm
car_brbo_norm<-brbo_norm
car_brbo_norm$forbin<-as.factor(car_brbo_norm$forbin)
levels(car_brbo_norm$forbin)<-c("transit", "forage")
car_brbo_norm$ID<-factor(car_brbo_norm$ID)
folds <- groupKFold(car_brbo_norm$ID) # should be same as non norm data

# resample method
train_control <- trainControl( method="LGOCV",
                               index=folds,
                               classProbs = TRUE,
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary)
# optimisation pars
mtry <- 4
tunegrid <- expand.grid(.mtry=mtry)

# raw run
rf_raw <- caret::train(forbin~sst_mn+sst_sd+chl_mn+chl_sd+mfront_sd+pfront_mn+ex_bathy+ex_slope,
                           data=car_brbo, method="rf", metric='ROC', 
                           tuneGrid=tunegrid, trControl=train_control)
Sys.time()
# norm run
rf_norm <- caret::train(forbin~sst_mn+sst_sd+chl_mn+chl_sd+mfront_sd+pfront_mn+ex_bathy+ex_slope,
                       data=car_brbo_norm, method="rf", metric='ROC', 
                       tuneGrid=tunegrid, trControl=train_control)
