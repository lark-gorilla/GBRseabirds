# Modelling (kernelUD approach) sent for running on Leeds HPC

library(dplyr)
library(ranger)
library(caret)
library(pROC)


set.seed(24) # so that random samples and random forest models can be reproduced

####~~ Tune hyperparameters of individual colony models and niave (all col) model ~~~####
sp_groups <- c('BRBO', 'MABO', 'RFBO', 'SOTE','WTST', 'WTLG','FRBD', 'TRBD', 'NODD', 'TERN')
for(k in sp_groups)
{
  for(j in 1:5)
  {  
  dat<-read.csv(paste0('kernhull_pts_sample/', k, '_kernhull_sample', j, '.csv'))
  dat$X<-NULL
  
  ## Tune individual colony models
    
  # hack to make inverse of groupKfolds()
  dat$temp<-1:nrow(dat)
  folds<-lapply(split(dat, dat$spcol), function(x){x$temp})
  
  train_control <- trainControl(method="LGOCV",index=folds,
                                classProbs = TRUE, savePredictions = TRUE,
                                summaryFunction = twoClassSummary, verboseIter = TRUE)
  
  tunegrid <- expand.grid(mtry=c(2:6),  splitrule = "gini", min.node.size = c(5,10,20,50))
  
  
  rf2 <- caret::train(x=dat[,c('sst','sst_sd','chl','chl_sd','mfr_sd', 'pfr_sd',
                               'mfr','pfr','bth','slp')],
                      y=dat[,'forbin'], method="ranger", num.trees=500, seed=24, metric='ROC', 
                      tuneGrid=tunegrid, trControl=train_control, verbose = TRUE)
  print(rf2)
  # add column to output predictions with colony for TEST partition
  rf2$pred$spcol<-unlist(lapply(folds, function(x){rep(dat[-x,]$spcol, nrow(tunegrid))}))
  
  indiv_aucz<-rf2$pred%>%group_by(Resample, spcol, mtry, min.node.size)%>%
    summarise(auc=as.double(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'), direction="<")$auc),
              thresh=coords(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
              sens=coords(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
              spec=coords(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
  indiv_aucz$TSS=indiv_aucz$sens+indiv_aucz$spec-1
  
  write.csv(indiv_aucz, paste0('col_tune/', k, '_indiv_col_tune', j, '.csv'), row.names=F, quote=F)
  rm(rf2)
  
  ## All colony niave model
 
  #folds2 <- groupKFold(dat$spcol)
  folds2<-lapply(folds, function(x){dat[-x,]$temp})# same as above but with colnames
  
  train_control <- trainControl( method="LGOCV",index=folds2,
                                 classProbs = TRUE, savePredictions = TRUE,
                                 summaryFunction = twoClassSummary, verboseIter = TRUE)
  
  rf_allcol <- caret::train(x=dat[,c('sst','sst_sd','chl','chl_sd','mfr_sd', 'pfr_sd',
                                      'mfr','pfr','bth','slp')],
                             y=dat[,'forbin'], method="ranger", num.trees=500, seed=24, metric='ROC', 
                             tuneGrid=tunegrid, trControl=train_control, verbose=T)
  # same tunegrid as indiv col models
  
  print(rf_allcol)

  allcol_aucz<-rf_allcol$pred%>%group_by(Resample, mtry, min.node.size)%>%
    summarise(auc=as.double(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'), direction="<")$auc),
              thresh=coords(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
              sens=coords(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
              spec=coords(pROC::roc(obs, Core, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
  allcol_aucz$TSS=allcol_aucz$sens+allcol_aucz$spec-1
  
  write.csv(allcol_aucz, paste0('col_tune/', k, '_all_col_tune', j,'.csv'), row.names=F, quote=F)
  rm(rf_allcol)
  }
}
