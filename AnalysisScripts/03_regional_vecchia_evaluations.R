# 03_regional_vecchia_evaluations.R
# Runs sampler on sub-regions of the data using the full Cholesky decomposition
# and the Vecchia approximation with m={10,25,50,100}. Code for producing
# Table 1 is found at the end of the script.
# IN: ../MCMC_Input/argo_data_january.RData
#     ../MCMC_Input/initial_parameters.RData
# OUT: Prints xtable object
# TEST: Computes predictions on coarser grid

test_mode=T
library(tidyverse)
library(BayesianOHC)

load('../MCMC_Input/argo_data_january.RData',verb=T)
load('../MCMC_Input/initial_parameters.RData',verb=T)

region_names=c("AtlN","PacNE","PacNW","AtlTrop","PacTropW","PacTropE",
               "Ind","PacSo","IndSo","AtlSo")
test_mode=T
if(test_mode){
  M=10
  global_grid=gen_masked_grid(10,10)
}else{
  M=1000
  global_grid=gen_masked_grid(1,1)
}
ms=c('chol',10,25,50,100)

yearlist=2007:2016
mean_field_values=array(NA,c(length(region_names),
                             length(ms),length(yearlist)))
trend_field_values=array(NA,c(length(region_names),
                              length(ms),length(yearlist)))
anomaly_field_values=array(NA,c(length(region_names),
                                length(ms),length(yearlist)))

myncores=5
#number of iterations to run:

for(regioni in 1:length(region_names)){
  region_name=region_names[regioni]
  region_data=in_region(argo_data_ordered,region_name)
  region_grid=in_region(global_grid,region_name)

  #needed for numeric integration
  scalarvec=gen_deltas_from_grid(region_grid)

  pred_df=region_data[,c(names(region_grid),'years')]
  pred_df$obs=T
  region_grid$obs=F
  for (year in 2007:2016) {
    region_grid$years=year
    pred_df=rbind(pred_df,region_grid)
  }
  for(mi in 1:length(ms)){
    m=ms[mi]
    print(paste('region_name:',region_name,'m:',m))

    if(m=="chol"){
      nyears=sapply(yearlist,function(year) sum(region_data$years==year))
      #trivial grouping (all obs grouped together) is equivalent to the full cholesky:
      grouping_list=lapply(nyears,groupinfo_full_chol)
      npred_years=sapply(yearlist,function(year) sum(pred_df$years==year))
      grouping_list_preddf=lapply(npred_years,groupinfo_full_chol)
    }else{
      m=as.numeric(m)
      grouping_list=compute_grouping_list(region_data,m=m,verb=F,
                                          ncores = myncores)
      grouping_list_preddf=compute_grouping_list(pred_df,m=m,verb=F,
                                                 ncores=myncores)
    }
    #supressing output is convenient if running in a loop:
    sink('/dev/null')
    region_sampler_out=run_mcmc_sampler(region_data,initparams,M=M,m=m,
                                        grouping_list=grouping_list,
                                        grouping_list_preddf = grouping_list_preddf,
                                        region_name=region_name,ncores = myncores)
    sink()
    posterior_densities=unlist(sapply(region_sampler_out,function(x) x$current_likelihood))+
      unlist(sapply(region_sampler_out,function(x) x$current_prior))
    map_params=region_sampler_out[[which(posterior_densities==max(posterior_densities))[1]]]$stored_params

    region_augdata=augment_data(region_data,map_params)
    krig_out=compute_lincomb_predictions(obsdata = region_augdata,
                                         myparams=map_params,
                                         pred_locs=region_grid,
                                         grouping_list_preddf=grouping_list_preddf,
                                         ret_var=F,scalarvec = scalarvec,
                                         ncores = myncores)
    mean_field_values[regioni,mi,]=krig_out$mu0
    trend_field_values[regioni,mi,]=krig_out$slope
    anomaly_field_values[regioni,mi,]=krig_out$resid
  }
}

###Table 1 Code

fractional_errors<-function(x){
  #sum over regions for global integral and take median over years
  return(x%>%apply(c(1,2),sum)%>%
           apply(1,function(x) abs(x/x[1]-1))%>%
           apply(1,median))
}

mean_tablevals=fractional_errors(mean_field_values)
trend_tablevals=fractional_errors(trend_field_values)
anom_tablevals=fractional_errors(anomaly_field_values)

table1_dat=data.frame(m=ms,
  "mean_field"=mean_tablevals,
  "trend_field"=trend_tablevals,
  "anom_field"=anom_tablevals)
require(xtable)
xtable(table1_dat,digits=-6)

