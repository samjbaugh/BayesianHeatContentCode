# 02_get_initial_parameters.R
# Calculates initial parameter fields by using a moving window approach where
# MLEs are computed in a window around each gridcell, and then hyperparameter
# surfaces are fit to the results. Code for Table 2 is at the end of the script.
# IN:  ../MCMC_Input/argo_data_january.RData
# OUT: ../MCMC_Input/knot_points.RData
#      ../MCMC_Input/initial_parameters.RData
# TEST: Does not save output, uses coarser grid for moving window calculations

testmode=T
library(tidyverse)
library(BayesianOHC)
library(pbmcapply)

load('../MCMC_Input/argo_data_january.RData',verb=T)
knot_points=gen_masked_grid(latres=8,lonres=16)
save(knot_points,file='../MCMC_Input/knot_points.RData')

verb=F

if(testmode){
  param_grid=gen_masked_grid(latres=24,lonres=24)
}else{
  param_grid=gen_masked_grid(latres=6,lonres=6)
}
ngrid=dim(param_grid)[1]

get_MLEs_window<-function(gridi)
{
  mydata=argo_data_ordered
  retval<-tryCatch(
    expr={
    mylon=param_grid$lon_degrees[gridi]
    mylat=param_grid$lat_degrees[gridi]
    window_size=10
    in_window=mygc_dist_degrees(mylon,mydata$lon_degrees)<window_size &
      abs(mydata$lat_degrees-mylat)<window_size
    in_window_data=mydata[in_window,]
    a=get_stationary_MLEs(in_window_data,verb=verb,est_nugget = T)},
    error=function(e) {print(e); return(NULL)})
    return(retval)
}

mycores=5
grid_output=pbmclapply(1:ngrid,get_MLEs_window,mc.cores=mycores)
good_vals=!sapply(grid_output,is.null)
grid_output=grid_output[good_vals]
moving_window_out_df=do.call(rbind,lapply(grid_output,data.frame))%>%
  cbind(param_grid[good_vals,])

inithypers=list()
initparams=list(basis_fields=knot_points,
                hyperparam_list=list(),
                basis_mu=knot_points)

for(varname in c('theta_lat','theta_lon','nugget','phi')){
  linkfun=default_linkfuns[[varname]]
  invlinkfun=default_invlinkfuns[[varname]]
  
  moving_window_out_df$z=invlinkfun(moving_window_out_df[[varname]])

  #estimate stationary hyperparameters
  esthypers=get_stationary_MLEs(mydata=moving_window_out_df %>%
                                  mutate(years=0),
                                yearlist=c(0),est_mu=T,est_slope = F,
                                est_nugget=F,iso=T)

  range=esthypers$theta_lat
  inithypers[[varname]]=list('range'=range,'mu'=esthypers$mu0,
                             'phi'=esthypers$phi,'nugget'=esthypers$nugget,
                             'theta_lat'=range,'theta_lon'=range)

  to_krig_df=moving_window_out_df %>% mutate(z=invlinkfun(moving_window_out_df[[varname]]),obs=T) %>%
    dplyr::select('lat_rad','lon_rad','lat_degrees','lon_degrees','z','obs') %>%
    rbind(knot_points %>% dplyr::select('lat_rad','lon_rad','lat_degrees','lon_degrees') %>%
            mutate('z'=NA,'obs'=F)) %>%
    mutate(theta_lon=esthypers$theta_lon,theta_lat=esthypers$theta_lat,
           phi=esthypers$phi,nugget=esthypers$nugget,mu=esthypers$mu0)
  preds=linkfun(krig_argo_field(to_krig_df,pred_indices=!to_krig_df$obs)$pred)

  to_krig_df$z[!to_krig_df$obs]=preds
  knot_params=to_krig_df%>%dplyr::filter(!obs)

  ####Convert to basis
  nknot=dim(knot_params)[1]
  c_knot_knot_mat=outer(1:nknot,1:nknot,cyl_cor_single,knot_params)
  knot_points[[paste('basis',varname,sep='_')]]=
    solve(t(chol(c_knot_knot_mat)),
          (invlinkfun(preds)-inithypers[[varname]]$mu)/
            sqrt(inithypers[[varname]]$phi))

  initparams$basis_fields=knot_points
  initparams$hyperparam_list[[varname]]=inithypers[[varname]]
}

#do mean/trend fields separately by computing joint marginal posterior means
varname='slope'
moving_window_out_df$z=(moving_window_out_df[[varname]])
esthypers_slope=get_stationary_MLEs(mydata=moving_window_out_df %>%
                                   mutate(years=0),
                                 yearlist=c(0),est_mu=F,est_slope = F,
                                 est_nugget=F,iso=T)
range=esthypers_slope$theta_lat
initparams$hyperparam_list[[varname]]=list('range'=range,'mu'=esthypers_slope$mu0,
                                           'phi'=esthypers_slope$phi,'nugget'=esthypers_slope$nugget,
                                           'theta_lat'=range,'theta_lon'=range)

varname='mu0'
moving_window_out_df$z=(moving_window_out_df[[varname]])
esthypers_mu=get_stationary_MLEs(mydata=moving_window_out_df %>%
                                   mutate(years=0),
                                 yearlist=c(0),est_mu=T,est_slope = F,
                                 est_nugget=F,iso=T)
initparams$hyperparam_list[[varname]]=
  list('range'=range,'mu'=esthypers_mu$mu0,
       'phi'=esthypers_mu$phi,'nugget'=esthypers_mu$nugget,
       'theta_lat'=range,'theta_lon'=range)

varname_list=names(linkfuns)
augdata_temp=augment_data(argo_data_ordered%>%mutate(mu0=0,slope=0),
                          initparams,
                         varname_list[!(varname_list%in%c('mu0','slope'))],
                         linkfuns)

grouping_list=compute_grouping_list(augdata_temp,m=50)
veccmat_list=build_veccmat_list_grouped(augdata_temp,grouping_list)

initparams=sample_mean_trend(augdata_temp,initparams,
                             veccmat_list,yearlist = 2007:2016)

augdata_init=augment_data(augdata_temp,initparams,
                         c('mu0','slope'),linkfuns)

if(!testmode){
  save(initparams,augdata_init,
       file='../MCMC_Input/initial_parameters.RData')
}

###Table 2 Code:

load('../MCMC_Input/initial_parameters.RData',verb=T)
k1=do.call(rbind,lapply(initparams$hyperparam_list[1:4],
                        function(x) data.frame(range=convert_theta_lat_to_effective_range_deg(x$range),
                                               q16=qlnorm(.25,x$mu,sqrt(x$phi)),
                                               mu_norm=exp(x$mu),
                                               q84=qlnorm(.75,x$mu,sqrt(x$phi)))))
k1['phi',-c(1)]=sqrt(k1['phi',-c(1)])
k1['theta_lat',-c(1)]=convert_theta_lat_to_effective_range_deg(k1['theta_lat',-c(1)])
k1['theta_lon',-c(1)]=convert_theta_lat_to_effective_range_deg(k1['theta_lon',-c(1)])

k2=do.call(rbind,lapply(initparams$hyperparam_list[5:6],function(x) data.frame(
  range=convert_theta_lat_to_effective_range_deg(x$range),
  q16=qnorm(.25,x$mu,sqrt(x$phi)),
  mu_norm=x$mu,
  q84=qnorm(.75,x$mu,sqrt(x$phi)))))
kk=rbind(k1,k2) %>% cbind(data.frame(varname=rownames(.)),.)
Table2=kk %>% mutate(Dist=c(rep("Log-Normal",4),rep("Normal",2))) %>%
  mutate(units=c("Degrees","Degrees","Unitless","GJ/m2","GJ/m2","GJ/(m2year)"))
Table2
