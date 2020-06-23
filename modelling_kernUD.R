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

#set.seed(7) really important to set.seed with RF??

# read in data
brbo<-read.csv('C:/seabirds/data/BRBO_forKern_modelready.csv')

hulls<-read.csv('C:/seabirds/data/BRBO_hulls_modelready.csv')

# clean up data

# Collapse ID to just remove owner and species, this combines all Soanes Dog
# and Swains cols - check for differences in interval
hulls$ID<-do.call(c, lapply(strsplit(as.character(hulls$ID), '_'), function(x)x[3]))
hulls$ID<-substr(hulls$ID, 1, (nchar(hulls$ID)-4)) 

# remove sitting behaviour and error pts

brbo$forbin=1
hulls$forbin=0

brbo<-rbind(brbo, hulls)


# set bathy vals >0 to NA then make positive
brbo[brbo$ex_bathy>0,]$ex_bathy<-NA
brbo$ex_bathy<-sqrt(brbo$ex_bathy^2)
brbo<-na.omit(brbo)

brbo%>%group_by(ID)%>%summarise(nP=length(which(forbin==1)),nA=length(which(forbin==0)))

brbo<-brbo%>%group_by(ID)%>%mutate(npres=length(which(forbin==1)))

brbo<-rbind(brbo%>%filter(forbin==1), brbo%>%filter(forbin==0)%>%group_by(ID)%>%
  sample_n(size=(unique(npres)*3), replace=T)) # 3:1 Pa to A

brbo$npres<-NULL

# niche overlap
enviro_std<-decostand(brbo[,c(2:11)], method="standardize")
enviro_rda<-rda(enviro_std, scale=T)
screeplot(enviro_rda)
enviro.sites.scores<-as.data.frame(scores(enviro_rda, choices=1:4, display='sites', scaling=1))
enviro.sites.scores<-data.frame(enviro.sites.scores,brbo[,c(1,12)])
ggplot(data=enviro.sites.scores, aes(x=PC1, y=PC2))+
  geom_density_2d(aes(colour=factor(forbin)))+facet_wrap(~ID)+theme_bw()


# Remember nearshore front values which ==0 should be NA

# scale: normalise (0-1) each covariate per dataID
# Remember to check for outliers prior to normaliziation
normalized<-function(x){(x-min(x))/(max(x)-min(x))}
brbo_norm<-brbo %>% group_by(ID) %>% mutate_if(is.numeric, normalized)%>%as.data.frame()

for( i in unique(brbo_norm$ID))
{
  #pr1<-table(brbo_norm[brbo_norm$ID==i,]$forbin)[1]/table(brbo_norm[brbo_norm$ID==i,]$forbin)[2]
  cw<-ifelse(brbo_norm[brbo_norm$ID==i,]$forbin==0, 1, 3)  
  
  rf1<-ranger(factor(forbin)~sst_mn+sst_sd+chl_mn+chl_sd+mfront_sd+pfront_mn+ex_bathy+ex_slope,
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
  
  # SPAC assessment
  sp1<-spline.correlog(brbo_norm[brbo_norm$ID==i,]$Longitude,
                       brbo_norm[brbo_norm$ID==i,]$latitude, residuals(rf1, type="pearson"),
                       na.rm=T,latlon=T,resamp=10)
  plot(sp1)
  
  
  
  print(i)
  print(Sys.time())
}

# check cols selected and size of matrix before running
br1<-as.data.frame(brbo_norm[,c(1, 12, 13:28)])%>%gather(variable, value, -ID, -forbin)
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


#caret approach

#test_dat<-as.data.frame(brbo_norm[brbo_norm$ID %in% unique(brbo_norm$ID)[c(2,4,5)],])

#try subsample data to something that will actually run
brbo_norm%>%group_by(ID)%>%summarise(n())
test_dat<-brbo_norm%>%group_by(ID)%>%
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
