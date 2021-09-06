# 14_compute_lowo_cv.R
# Computes windowed cross validation (also leave-one-window-out or LOWO) on 
# each of the models model names are 'fullmomdel', 'stat_phi', 'stat_theta_lat', 
# 'stat_theta_lon', 'stat_nugget', 'stat_all', 'chordal' (GpGp), 
# 'spacetime' (GpGp), and 'chordal_nonstat'
# IN: Reads the output for each model; see script for exact locations
# OUT: For each model as named above saves an object of the form:
#      ../ValidationData/lowodf_<model_name>.RData
#      What is saved is a dataframe with fields "validval" and "validsd" for
#      the cross-validation predicted values and standard deviations for
#      each location
# TEST: N/A, be careful this script will overwrite existing data

library(tidyverse)
library(BayesianOHC)

cv_type='lowo'
cv_options=list(lowo_width=2)
myncores=10

####Fullmodel CV
{
  model_name='fullmodel'
  load('../MCMC_Output/Output_20k_run/map_configuration.RData',verb=T)
  lowodf_fullmodel=compute_cv_df(map_augdata,cv_params=veccmat_list,ncores=myncores,
                                  cv_type=cv_type,cv_options=cv_options)
  save(lowodf_fullmodel,file=paste('../ValidationData/lowodf_',model_name,'.RData',sep=''))
}

####Statmodel CV
{
  statmodel_names=c('stat_phi','stat_theta_lat','stat_theta_lon',
                    'stat_nugget','stat_all')
  for(model_name in statmodel_names){
    print(model_name)
    load(paste('../ValidationData/map_params_',model_name,'.RData',sep=''),verb=T)
    lowodf_statmod=compute_cv_df(augdata_statmod,cv_params=veccmat_list_statmod,ncores=myncores,
                                  cv_type=cv_type,cv_options=cv_options)
    save(lowodf_statmod,file=paste('../ValidationData/lowodf_',model_name,'.RData',sep=''))
  }
}

####GpGp Fits
{
  for(model_name in c('spacetime','chordal')){
    load(file=paste('../ValidationData/gpgp_',model_name,'_cvparams.RData',sep=''),verb=T)
    lowo_df_gpgp=compute_cv_df(augdata_gpgp,
                               cv_params=cv_params,ncores=myncores,
                               cv_type=cv_type,cv_options=cv_options)
    save(lowo_df_gpgp,file=paste('../ValidationData/lowodf_',model_name,'.RData',sep=''))
  }
}

####NonStat Chordal
{
  model_name='chordal_nonstat'
  load('../ValidationData/map_sample_chordal_nonstat.RData',verb=T)
  lowodf_chordal_nonstat=compute_cv_df(augdata_ordered,cv_params=veccmat_list,ncores=myncores,
                                        cv_type=cv_type,cv_options=cv_options)
  save(lowodf_chordal_nonstat,file=paste('../ValidationData/lowodf_',model_name,'.RData',sep=''))

}


####Levitus
{
  load('../MCMC_Input/argo_data_january.RData')
  load('../ValidationData/levitus_anomalies.RData',verb=T)
  yearlist=2007:2016
  require('geodist')
  distmats=lapply(yearlist,function(year)
    geodist(subset(argo_data_ordered,years==year)[,c('lon_degrees','lat_degrees')],
            measure='haversine')/1e3)
  lowodf_levitus=compute_cv_df(argo_data_levitus_anoms,
                                cv_function=levitus_preds,yearlist=yearlist,
                                cv_params=distmats,ncores=1,
                                cv_type=cv_type,cv_options=cv_options)
  save(lowodf_levitus,file='../ValidationData/lowodf_levitus.RData')
}


###NSGP (not run for paper)
if(F){
  load('../MCMC_Input/argo_data_january.RData')
  load('../ValidationData/nspg_samples.RData',verb=T)
  load('../ValidationData/knot_points_chordal.RData',verb=T)

  burnin=5000
  M=dim(samples)[1]
  nsgp_post_means=data.frame(t(apply(samples[burnin:M,],2,mean)))
  argo_data_chord=argo_data_ordered%>%add_3d_coord

  load(file='../ValidationData/nsgp_grouping_list')

  out=convert_nsgp_to_cvparams(nsgp_post_means,knot_points,argo_data_chord,
                               grouping_list)
  augdata_nsgp=out$augdata
  cv_params=out$cv_params

  lowodf_nsgp_fullnonstat=compute_cv_df(augdata_nsgp,cv_params=cv_params)
  save(lowodf_nsgp_fullnonstat,file='../ValidationData/lowodf_nsgp_fullnonstat.RData')
}


