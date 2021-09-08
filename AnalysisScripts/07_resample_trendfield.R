# 07_resample_trendfield.R
# Obtains samples of the posterior distribution of the trend field by
# resampling from the conditional posterior distributions at each iteration
# IN: ../MCMC_Input/argo_data_january.RData
#     ../MCMC_Output/Output_20k_run/sample_matrices.RData
#     ../MCMC_Output/Output_20k_run/map_configuration.RData
# OUT: ../MCMC_Output/Output_20k_run/trend_resampled.RData
#     ../MCMC_Output/Output_20k_run/trend_resamp_intonly.RData
# TEST: Does not overwrite data

test_mode=T
library(dplyr)
library(BayesianOHC)

load('../MCMC_Input/argo_data_january.RData',verb=T)
load('../MCMC_Output/Output_20k_run/sample_matrices.RData',verb=T)
load('../MCMC_Output/Output_20k_run/map_configuration.RData',verb=T)
mygrid=gen_masked_grid(1,1)
#replace with samples from each iteration:
resampled_slope_list=list()
nresamp=100

resample_slope<-function(iter){
  myhypermu=sample_matrices$hyper_mu0_mean[iter]
  mymu0=sample_matrices$mu0[,iter]
  myslope=sample_matrices$slope[,iter]
  mu_sampleout$retparams$hyperparam_list$mu0$mu=myhypermu
  mu_sampleout$retparams$basis_mu$basis_mu0=mymu0
  mu_sampleout$retparams$basis_mu$basis_slope=myslope
  mu_sampleout$mean=c(myhypermu,myslope,mymu0)
  return(resample_slopes(mu_sampleout,nresamp=nresamp,mygrid))
}

load(file='../MCMC_Output/Output_20k_run/mu_slope_df.RData',verb=T)
nresamp=100
#sample integrated trend values, thinning by 10
thinned_iters=seq(1,dim(mu_slope_df)[1],by=10)
resamp_int_trends=sapply(thinned_iters,
                         function(i) rnorm(nresamp,
                                           mu_slope_df$slope[i],
                                           mu_slope_df$slope_sd[i]))
slope_ci=quantile(resamp_int_trends,c(.05,.5,.95))

#coarser thin, only need to find close match to integrated value:
thinned_iters_coarse=seq(1,dim(mu_slope_df)[1],by=1000)

resampled_slope_list=parallel::mclapply(thinned_iters_coarse,
                                        resample_slope,mc.cores=5)
all_field_samples=do.call(rbind,lapply(resampled_slope_list,function(x) x$trend_samples))
all_int_means=do.call(c,lapply(resampled_slope_list,function(x) x$intslope_samples))

up_quant_i=which(abs(all_int_means-slope_ci[1])==
            min(abs(all_int_means-slope_ci[1])))
low_quant_i=which(abs(all_int_means-slope_ci[3])==
            min(abs(all_int_means-slope_ci[3])))

if(!test_mode){
  save(resamp_int_trends,all_field_samples,all_int_means,up_quant_i,low_quant_i,
       file='../MCMC_Output/Output_20k_run/trend_resampled.RData')
  save(resamp_int_trends,file='../MCMC_Output/Output_20k_run/trend_resamp_intonly.RData')
}



