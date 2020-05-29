# Modelling code

library(dplyr)
library(ggplot2)
library(GGally)
library(randomForest)
library(caret)

#set.seed(7) really important to set.seed with RF??

# read in data
brbo<-read.csv('C:/seabirds/data/tracking_modelready.csv')

hulls<-read.csv('C:/seabirds/data/BRBO_hulls_modelready.csv')


# clean up data

# remove sitting behaviour and error pts

brbo<-filter(brbo, embc %in% c('commuting', 'relocating', 'foraging') )
brbo$forbin<-1
brbo[brbo$embc=='commuting',]$forbin<-0

# remove NA pts (from chl, sst and front land mask mainly )
brbo$trip_id<-NULL # temporary fix for NICH trip_ID==NA
brbo<-na.omit(brbo)

# Remember nearshore front values which ==0 should be NA
# Remember some bad tracks that should be filtered prior to modelling
# rememeber could remove bathy>0 vals

# scale: normalise (0-1) each covariate per dataID
# Remember to check for outliers prior to normaliziation
normalized<-function(x){(x-min(x))/(max(x)-min(x))}
brbo_norm<-brbo %>% group_by(ID) %>% mutate_if(is.numeric, normalized)
brbo_norm$forbin<-brbo$forbin
# check correlation
#ggpairs(brbo_norm[,3:14]) # takes ages

cor(brbo_norm[,4:14])
# chl and bathy cor = 0.54 but keep
# fronts metrics cor. Keep pfront_mn and mfront_sd @ 0.47

for( i in unique(brbo_norm$ID))
{
rf1<-randomForest(factor(forbin)~sst_mn+sst_sd+chl_mn+chl_sd+mfront_sd+pfront_mn+ex_bathy+ex_slope,
                  data=brbo_norm[brbo_norm$ID==i,], ntree=500, mtry=4, importance=T)
print(rf1)
print(varImpPlot(rf1))
brbo_norm$p1<-predict(rf1, newdata=brbo_norm, type='prob')[,2] # prob of foraging 0-1
names(brbo_norm)[which(names(brbo_norm)=='p1')]<-substr(i, 1, (nchar(i)-4))
print(i)
print(Sys.time()) 
}

#caret approach

test_dat<-as.data.frame(brbo_norm[brbo_norm$ID %in% unique(brbo_norm$ID)[c(2,4,5)],])

# make forbin a factor with labels attached rather than 0, 1
test_dat$forbin<-as.factor(test_dat$forbin)
levels(test_dat$forbin)<-c("transit", "forage")

# do somne variable vis
ggplot(data=test_dat, aes(x=sst_mn, colour=forbin))+
  geom_density()

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
mtry <- 1:4
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- caret::train(forbin~sst_mn+sst_sd+chl_mn+chl_sd+mfront_sd+pfront_mn+ex_bathy+ex_slope,
                    data=test_dat, method="rf", metric='ROC', 
                    tuneGrid=tunegrid, trControl=train_control)
print(rf_default)

rf_default$resample

head(rf_default$pred)

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


 


# try caret approach
control <- trainControl(method="repeatedcv", number=10, repeats=3, indicies=folds)
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(x))
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(Class~., data=dataset, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)
