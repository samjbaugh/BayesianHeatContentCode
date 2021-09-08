# 04_run_global_mcmc.R
# Runs MCMC sampler for non-stationary cylindrical distance model
# IN: ../MCMC_Input/argo_data_january.RData
# OUT: Writes output to directory
#       '../MCMC_Output/Output_run_label' where run_label can be customized;
#       samples are stored in the "Samples" subdirectory and at the end MAP
#       parameters/fields are determined and stored. Besides samples, includes 
#       the following output:
#             ../MCMC_Output/Output_run_label/map_configuration.RData
#             ../MCMC_Output/Output_run_label/ohc_df.RData and 
#             ../MCMC_Output/Output_run_label/mu_slope_df.RData 
#       in the output directory
# TEST: Runs sampler for a small number of iterations on a subset of the data/
#       running in test mode will write output to directory suffixed with 
#       "test_run".

test_mode=T
library(tidyverse)
library(BayesianOHC)

myncores=5
m=50

if(test_mode){
  M=200
  load('../MCMC_Input/argodata_subset.RData',verb=T)
  mcmc_input_data=argo_data_subset
  ohc_grid=gen_masked_grid(24,24)
  run_label=paste('testrun_',substr(date(),5,7),
                  substr(date(),9,10),'_',substr(date(),21,24),sep='')
  save_iter=10
}else{
  M=20000
  load('../MCMC_Input/argo_data_january.RData',verb=T)
  mcmc_input_data=argo_data_ordered
  ohc_grid=gen_masked_grid(2,2)
  #previous run with run_label=20k_run
  run_label=paste('20k_rerun_',substr(date(),5,7),
                  substr(date(),9,10),'_',substr(date(),21,24),sep='')
  save_iter=100
}
load('../MCMC_Input/initial_parameters.RData',verb=T)

#Setting save_iter means that
region_sampler_out=run_mcmc_sampler(ordered_data=mcmc_input_data,
                                    initparams=initparams,
                                    M=M,m=m,pred_locs=ohc_grid,
                                    save_iter=save_iter,pred_iter=100,
                                    run_label=run_label,
                                    outdir='../MCMC_Output',
                                    ncores = myncores)

#change label if want to use new run
run_label='Output_20k_run'
outdir=paste('../MCMC_Output/',run_label,sep='')
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
  save(map_augdata,map_sample,grouping_list,veccmat_list,mu_sampleout,
       file=paste('../MCMC_Output/',run_label,'/map_configuration.RData',sep=''))
}

