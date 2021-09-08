# 10_GpGp_stationary_chordal_fit.R
# Fits fully stationary chordal distance model using the GpGp package. Also
# converts fitted models into "cv_params" object that can be passed into
# cross validation functions
# IN: ../MCMC_Input/argo_data_january.RData
# OUT: ../ValidationData/gpgp_chordal_modelfits.RData
#      ../ValidationData/gpgp_chordal_cvparams.RData
# TEST: Does not overwrite output

test_mode=T
library(tidyverse)
library(BayesianOHC)

load('../MCMC_Input/argo_data_january.RData',verb=T)
yearlist=2007:2016

###Fit models:
model_fit_list=list()
for(myyear in yearlist){
  yeardata_ordered=subset(argo_data_ordered,years==myyear)
  X=model.matrix(~lat_degrees*lon_degrees+
                   I(lat_degrees^2)+I(lon_degrees^2),
                 data=yeardata_ordered)
  locs_lonlat=as.matrix(yeardata_ordered[,c('lon_degrees','lat_degrees')])

  NNarray=GpGp::find_ordered_nn(locs_lonlat,m=50,lonlat=T)
  gpgp_model_fit=GpGp::fit_model(y=yeardata_ordered$vhc_obs,
                           X=X,
                           locs=locs_lonlat,m_seq = c(10,50),
                           NNarray = NNarray,
                           covfun_name="matern_sphere")
  gpgp_model_fit$NNarray=NNarray
  model_fit_list[[myyear-2007+1]]=gpgp_model_fit
}

if(!test_mode){
  save(model_fit_list,file='../ValidationData/gpgp_chordal_modelfits.RData')
}

###Compute parameters for cross validation functions:
load(file='../ValidationData/gpgp_chordal_modelfits.RData',verb=T)

augdata_gpgp=map_dfr(yearlist,~add_means(subset(argo_data_ordered,years==.),
                                         model_fit_list[[.-2006]]))

matern_wrapper_sphere<-function(covparms,df){
  covparms[1]=1
  GpGp::matern_sphere(covparms,as.matrix(df[,c('lon_degrees','lat_degrees')]))
}

#convert to cv_params object, needs to be done in order to make the model
#output compatible with the cross validation functions
cv_params=convert_gpgp_fit_to_cvparams(model_fit_list,augdata_gpgp,
                                       corrfun=matern_wrapper_sphere,
                                       yearlist=yearlist)

if(!test_mode){
  save(augdata_gpgp,cv_params,file='../ValidationData/gpgp_chordal_cvparams.RData')
}
