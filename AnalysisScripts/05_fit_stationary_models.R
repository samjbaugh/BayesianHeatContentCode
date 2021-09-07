# 05_fit_stationary_models.R
# Runs the MCMC sampler to fit models stationary in each of the parameter fields
# for the cross validation study
# IN: '../MCMC_Output/Output_20k_run/map_configuration.RData' (this is for
#     initializing the sampler; can also be initialized from the last iteration
#     of the full sampler or from some other configuration)
# OUT: For each of the models fit, writes files of the name
#      map_params_<modelname>.RData in the ValidationData sub-directory
#      containing the MAP parameters, fields, and veccmat objects
# TEST: Runs the sampler for a fewer number of iterations on subsetted data

testmode=T
library(tidyverse)
library(BayesianOHC)

#use map configuration from full run for initialization:
load(file=paste('../MCMC_Output/Output_20k_run/map_configuration.RData',sep=''),verb=T)
myncores=5
m=50

if(test_mode){
  M=10
  load('../MCMC_Input/argodata_subset.RData',verb=T)
  ordered_data=argo_data_subset
}else{
  M=5000
  load('../MCMC_Input/argo_data_january.RData',verb=T)
  ordered_data=argo_data_ordered
}

model_list=list('stat_phi'=list('stationary_params'=c('phi'),'iso'=F),
                'stat_theta_lat'=list('stationary_params'=c('theta_lat'),'iso'=F),
                'stat_theta_lon'=list('stationary_params'=c('theta_lon'),'iso'=F),
                'stat_nugget'=list('stationary_params'=c('nugget'),'iso'=F),
                'stat_all'=list('stationary_params'=c('phi','theta_lat',
                                         'theta_lon','nugget'),'iso'=F))
model_names=names(model_list)

fullmodel_params=map_sample$stored_params
map_augdata=augment_data(pred_locs=ordered_data,
                         myparams=map_sample$stored_params)

grouping_list=compute_grouping_list(ordered_data,m=m,verb=F,
                                    ncores = myncores)
for(model_name in model_names){
  stationary_initparams=convert_params_to_statiso(fullmodel_params,map_augdata,
                                                  model_list[[model_name]])

  print(paste('Fitting model:',model_name))
  region_sampler_out=run_mcmc_sampler(ordered_data,stationary_initparams,
                                      M=M,m=m,grouping_list=grouping_list,
                                      ncores=myncores)

  posterior_densities=
    unlist(sapply(region_sampler_out,function(x) x$current_likelihood))+
    unlist(sapply(region_sampler_out,function(x) x$current_prior))
  nonstat_map_params=region_sampler_out[[
    which(posterior_densities==max(posterior_densities))[1]]]$stored_params

  augdata_statmod=augment_data(ordered_data,nonstat_map_params)
  veccmat_list_statmod=build_veccmat_list_grouped(augdata_statmod,
                                                  grouping_list,ncores=5)

  if(!test_mode){
      save(nonstat_map_params,augdata_statmod,veccmat_list_statmod,
       file=paste('../ValidationData/map_params_',model_name,'.RData',sep=''))
  }
}
