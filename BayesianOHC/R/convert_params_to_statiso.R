#' Convert nonstationary anisotropic parameters to stationary or isotropic
#'
#' @param myparams Parameters to convert
#' @param myaugdata Augmented data to use in the conversion
#' @param mymodel Contains list of parameters to convert to stationarity and a
#' boolean "iso" for whether or not the model should be converted to isotropic
#' @export
convert_params_to_statiso<-function(myparams,myaugdata,mymodel){
  myparams$stationary_params=list()
  stationary_varnames=mymodel$stationary_params
  iso=mymodel$iso
  if(iso){
    if('theta' %in% stationary_varnames){
      myparams$stationary_params[['theta']]=median(myaugdata$theta_lat)
    }else{
      myparams$basis_fields$basis_theta=myparams$basis_fields$basis_theta_lat
    }
    myparams$hyperparam_list[['theta']]=myparams$hyperparam_list$theta_lat
    myparams$hyperparam_list$theta_lat=NULL
    myparams$hyperparam_list$theta_lon=NULL
    myparams$basis_fields$basis_theta_lat=NULL
    myparams$basis_fields$basis_theta_lon=NULL
  }
  for(varname in stationary_varnames){
    if(varname!='theta'){
      basis_varname=paste('basis',varname,sep='_')
      myparams$stationary_params[[varname]]=
        median(myaugdata[[varname]])
      myparams$basis_fields[[basis_varname]]=NULL
    }
  }
  return(myparams)
}
