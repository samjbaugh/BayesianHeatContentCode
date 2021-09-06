# 01_preprocess_argo_data.R
# Reads in the raw Argo data file (which contains the temperatures,
# pressures, and salinities  at each depth level) and computes numerical
# integrals with  respect to depth to calculate the vertical ocean heat
# content values at each latitude/longitude/time observation.
# IN: RawData/argo_object_01.RData
# OUT: MCMC_Input/argo_data_january.RData
# TEST: Does not save output

test_mode=T
library(BayesianOHC)
library(tidyverse)
library(pbapply)

#'01' is for january:
month_id=01
sw_sh<<-3989.411
sw_density<<-1028.319

raw_data_filename=paste('../RawData/argo_object_',
                        formatC(month_id,width=2,flag="0"),'.RData',sep='')
if(file.exists(raw_data_filename))
{
  load(raw_data_filename,verb=T)
}else{
  print('No file')
}

max_depth=unlist(lapply(mat_object$profPresAggr,function(x) max(x[[1]])))
lon_unshifted=mat_object$profLongAggr[1,]
lon_shifted=ifelse(lon_unshifted>180,lon_unshifted-360,lon_unshifted)

profile_df=data.frame(lat_degrees=mat_object$profLatAggr[1,],lon_degrees=lon_shifted,
                      lat_rad=deg2rad(mat_object$profLatAggr[1,]),lon_rad=deg2rad(lon_shifted),
                      max_depth=max_depth,years=mat_object$years,
                      days=mat_object$days,months=mat_object$months,hours=mat_object$hours,
                      minutes=mat_object$minutes,seconds=mat_object$seconds,
                      time=mat_object$days*24*60*60+
                        mat_object$hours*60*60+
                        mat_object$minutes*60+
                        mat_object$seconds,
                      float_id=mat_object$profFloatIDAggr[1,])

mm=length(mat_object$days)

calculate_vhc<-function(ii)
{
  obs_pressures=mat_object$profPresAggr[ii][[1]][[1]]
  obs_pressures_clean=obs_pressures[obs_pressures>0 & obs_pressures<2000]
  if(length(obs_pressures_clean)==0)
  {
    return(data.frame('vhc_obs'=NA))
  }
  obs_pressures_endpoints=c(0,obs_pressures_clean,2000)
  obs_temps=mat_object$profTempAggr[ii][[1]][[1]][obs_pressures>0 & obs_pressures<2000]
  obs_temps_endpoints=c(obs_temps[1],obs_temps,obs_temps[length(obs_temps)])
  vhc_obs=trapz(obs_pressures_endpoints,
                obs_temps_endpoints*sw_sh*sw_density)
  return(data.frame('vhc_obs'=vhc_obs))
}
vhc_df=do.call(rbind,pblapply(1:mm,calculate_vhc))
profile_df=cbind(profile_df,vhc_df)

#scale vhc_obs to giga-joules and time to days since Jan 1
argo_data_unordered=profile_df%>%
  dplyr::filter(max_depth>1900)%>%
  mutate(vhc_obs=vhc_obs/10^9)%>%
  mutate(time=(time-min(time))/(60^2*24))%>%
  mutate(z=vhc_obs)%>%
  dplyr::filter(.,!duplicated(.[,c('lat_degrees','lon_degrees','float_id')]))

argo_data_ordered=order_yeardata(argo_data_unordered)

if(!test_mode){
  save(argo_data_ordered,file='../MCMC_Input/argo_data_january.RData')
}
