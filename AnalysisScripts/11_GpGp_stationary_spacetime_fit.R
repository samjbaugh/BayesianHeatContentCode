# 11_GpGp_spacetime_fit.R
# Fits fully stationary spacetime model using the GpGp package. Also
# converts fitted models into "cv_params" object that can be passed into
# cross validation functions
# IN: ../MCMC_Input/argo_data_january.RData
# OUT: ../ValidationData/gpgp_spacetime_modelfits.RData
#      ../ValidationData/gpgp_spacetime_cvparams.RData
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
  locs_lonlattime=as.matrix(yeardata_ordered[,c('lon_degrees','lat_degrees',
                                                'time')])

  #spatio-temporal scaling determined through iterative process:
  st_scale=c(4.499155e-01,2.220130e+03)
  NNarray=GpGp::find_ordered_nn(locs_lonlattime,m=50,lonlat=T,st_scale=st_scale)
  gpgp_model_fit=GpGp::fit_model(y=yeardata_ordered$vhc_obs,
                                 X=X,
                                 locs=locs_lonlattime,m_seq = c(10,50),
                                 NNarray = NNarray,
                                 covfun_name="matern_spheretime")
  gpgp_model_fit$NNarray=NNarray
  model_fit_list[[myyear-2007+1]]=gpgp_model_fit
}

if(!test_mode){
  save(model_fit_list,file='../ValidationData/gpgp_spacetime_modelfits.RData')
}

###Compute parameters for cross validation functions:
load(file='../ValidationData/gpgp_spacetime_modelfits.RData',verb=T)

augdata_gpgp=map_dfr(yearlist,~add_means(subset(argo_data_ordered,years==.),
                                         model_fit_list[[.-2006]]))

matern_wrapper_spheretime<-function(covparms,df){
  covparms[1]=1
  GpGp::matern_spheretime(covparms,as.matrix(df[,c('lon_degrees','lat_degrees','time')]))
}

cv_params=convert_gpgp_fit_to_cvparams(model_fit_list,augdata_gpgp,
                                       corrfun=matern_wrapper_spheretime,
                                       yearlist=yearlist)

if(!test_mode){
  save(augdata_gpgp,cv_params,file='../ValidationData/gpgp_spacetime_cvparams.RData')
}
