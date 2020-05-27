# Modelling code

library(dplyr)
library(ggplot2)

# read in data
brbo<-read.csv('C:/seabirds/data/tracking_modelready.csv')

hulls<-read.csv('C:/seabirds/data/BRBO_hulls_modelready.csv')


# clean up data

# remove sitting behaviour and error pts

brbo<-filter(brbo, embc %in% c('commuting', 'relocating', 'foraging') )
brbo$forbin<-1
brbo[brbo$embc=='commuting',]$forbin<-0

# remove NA pts (from chl, sst and front land mask mainly )
brbo<-na.omit(brbo)

# Remember nearshore front values which ==0 should be NA
# Remember some bad tracks that should be filtered prior to modelling
# rememeber could reomve bathy>0 vals

# scale: normalise (0-1) each covariate per dataID
# Remember to check for outliers prior to normaliziation
normalized<-function(x){(x-min(x))/(max(x)-min(x))}
brbo_norm<-brbo %>% group_by(ID) %>% mutate_if(is.numeric, normalized)
