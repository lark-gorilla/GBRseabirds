# Random Forest stratified sampling playground
#https://github.com/imbs-hl/ranger/issues/142
#https://github.com/imbs-hl/ranger/issues/515
# check if adding case weights improves(https://github.com/imbs-hl/ranger/issues/167)

library(dplyr)
library(ranger)
library(pROC)
suppressWarnings()

samp_tree_func<-function(my.spcol='Rudi', my.dat=dat, colrep=3){
  ib1<-rep(0, nrow(my.dat))
  colsize<-length(which(my.spcol==my.dat$spcol))
  s1<-c(1:colsize, sample(1:colsize, size=colsize*colrep, replace=T))
  ib1[which(my.spcol==my.dat$spcol)]<-as.data.frame(table(s1))$Freq
  return(ib1)}

samp_tree_func_spread<-function(my.spcol='Rudi', my.dat=dat){
  ib1<-rep(0, nrow(my.dat))
  ib1[sample(1:nrow(my.dat), (nrow(my.dat)/2))]<-1
  colsize<-length(which(my.spcol==my.dat$spcol))
  s1<-c(1:colsize, sample(1:colsize, size=(nrow(my.dat)/2), replace=T))
  ib1[which(my.spcol==my.dat$spcol)]<-as.data.frame(table(s1))$Freq
  return(ib1)}

brbo_ntree<-read.csv('brbo_trans_ntree.csv')

sp_groups <- c('BRBO')
for(m in sp_groups)
{
  dat<-read.csv(paste0('kernhull_pts_sample/', m, '_kernhull_sample.csv'))
  dat$X<-NULL

for(h in unique(dat$spcol))
{
  tempdat<-dat[dat$spcol!=h,] # subset data to train colonies
  
  ncol<-length(unique(tempdat$spcol))
  
  temp_ntree<-sum(brbo_ntree[brbo_ntree$spcol%in%unique(tempdat$spcol),]$ntree)
  
  # create inbag: limit sampling (with replacement) of values for each tree to within each colony 
  inbag_list<-NULL
  for(k in unique(tempdat$spcol))
  {
    r1<- replicate(brbo_ntree[brbo_ntree$spcol==k,]$ntree, 
                   samp_tree_func_spread(my.spcol = k, my.dat=tempdat), simplify = FALSE)
    inbag_list<-c(inbag_list, r1)
  }
  
  #fit rfs
  rf_strat<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp ,
                   data=tempdat, num.trees=temp_ntree, inbag = inbag_list,
                   splitrule = "gini",  importance='impurity',probability =T, seed=24)
  
  #rf_simp<-ranger(forbin~sst+sst_sd+chl+chl_sd+mfr_sd+pfr_sd+pfr+mfr+bth+slp ,
  #                data=tempdat, num.trees=(ncol*100), 
  #                splitrule = "gini",  importance='impurity',probability =T, seed=24)
  
  # predict to dat
  dat$p1<-predict(rf_strat, data=dat)$predictions[,1]
  #dat$p2<-predict(rf_simp, data=dat)$predictions[,1]
  names(dat)[which(names(dat)=='p1')]<-paste0('stratWO_',h)
  #names(dat)[which(names(dat)=='p2')]<-paste0('simplWO_',h)
  print(h)
}

indiv_aucz<-tidyr::gather(dat[,c(1,2,16:ncol(dat))], Resample, pred, -forbin, -spcol)%>%
  group_by(Resample, spcol)%>%
  summarise(auc=as.double(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'), direction="<")$auc),
            thresh=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$threshold[1],
            sens=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$sensitivity[1],
            spec=coords(pROC::roc(forbin, pred, levels=c('PsuedoA', 'Core'),direction="<"),'best', best.method='youden', transpose=F)$specificity[1])
indiv_aucz$TSS=indiv_aucz$sens+indiv_aucz$spec-1

indiv_aucz$sp<-m

aucz_out<-rbind(aucz_out, indiv_aucz)
print(m)
}
write.csv(aucz_out, 'auc_strat_simpl_test2.csv', quote=F, row.names = F)
print(m)
