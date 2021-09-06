# 09_BayesianOHC_nonstat_chordal_fit.R
# Fits fully nonstationary (but isotropic) chordal fit using BayesianOHC
# IN: ../MCMC_Input/argo_data_january.RData
#     ../MCMC_Input/grouping_list.RData
#     ../MCMC_Output/Output_20k_run/map_configuration.RData
# OUT: Run_sampler writes to directory ../MCMC_Output/Output_chordal_nonstat
#      At the end MAP parameters for the cross validation study are stored in:
#      ../ValidationData/map_sample_chordal_nonstat.RData
# TEST: Does not overwrite data

testmode=T
library(tidyverse)
library(BayesianOHC)

myncores=5
m=50

load('../MCMC_Input/argo_data_january.RData',verb=T)
load('../MCMC_Input/grouping_list.RData',verb=T)
load('../MCMC_Output/Output_20k_run/map_configuration.RData',verb=T)

mcmc_input_data=argo_data_ordered
grouping_list2=compute_grouping_list(mcmc_input_data,m=50,ncores=myncores)

#start at the MAP values of the other parameters to minimize burn-in time:
map_params=map_sample$stored_params

#convert to isotropy:
chordal_initparams=convert_params_to_statiso(map_params,map_augdata,
                                             mymodel=list('stationary_params'=c(),
                                                          'iso'=T))

#convert units for range parameter, not sure how much this makes sense 
#but shouldn't matter after this parameter is re-burned in
ld2=function(x) {log(2*sin(x/2))}
chordal_initparams$hyperparam_list$theta$mu=
  map_params$hyperparam_list$theta_lat$mu%>%exp%>%
  convert_theta_lat_to_effective_range_deg%>%deg2rad%>%
  ld2

augdata_chordal=argo_data_ordered%>%
  augment_data_parallel(chordal_initparams,ncores=5)%>%
  add_3d_coord

if(test_mode){
  run_label='test_chordal_nonstat'
}else{
  run_label='chordal_nonstat'
}
outdir='~/MCMC_Storage/MCMC_Output'
region_sampler_out=run_mcmc_sampler(ordered_data=augdata_chordal,
                                    corrfun=chord_cor_single,
                                    initparams=chordal_initparams,
                                    grouping_list=grouping_list,
                                    M=10000,m=50,save_iter=100,
                                    run_label=run_label,
                                    outdir=outdir,
                                    ncores = myncores)

#change label if want to use new run
outdir=paste(outdir,run_label,sep='')
file_list=list.files(outdir)
max_posterior_density=-Inf
for(filename in file_list){
  load(paste(outdir,filename,sep=''),verb=T)
  posterior_densities=unlist(sapply(param_stores,function(x) x$current_likelihood))+
    unlist(sapply(param_stores,function(x) x$current_prior))
  filemax=max(posterior_densities)
  if(filemax>max_posterior_density){
    max_posterior_density=filemax
    map_sample=param_stores[[which(posterior_densities==filemax)[1]]]
  }
}

#save additional info for MAP (maximum a-posteriori) and last configuration:
map_augdata=augment_data(argo_data_ordered,map_sample$stored_params)
grouping_list=compute_grouping_list(argo_data_ordered,m=50)
veccmat_list=build_veccmat_list_grouped(map_augdata,
                                        grouping_list,ncores=1)
mu_sampleout=sample_mean_trend(map_augdata,map_sample$stored_params,veccmat_list,
                               return_full_posterior=T,ncores=1)

if(!test_mode){
  save(map_augdata,map_sample,veccmat_list,
       file='../ValidationData/map_sample_chordal_nonstat.RData')
}

