#' Augments data with parameter values
#' @description Adds parameters to a given data frame by kriging the basis values specified in
#' "myparams" to the locations in "mydata"
#'
#' @param pred_locs Locations for kriging the parameter fields
#' @param myparams Basis parameter values for kriging
#' @param varname_list List of parameter fields to krig
#' @param linkfuns List that contains an entry for each parameter field
#' specified in varname_list
#' @return Returns a data frame that is a copy of the input data but with new
#' parameter values for each of the fields specified in varname_list
#' @export
augment_data=function(pred_locs,myparams,
                      varname_list=NULL,
                      linkfuns=default_linkfuns){
  if(is.null(varname_list)){
    varname_list=names(myparams$hyperparam_list)
  }
  ret_augdata=pred_locs
  stationary_varnames=names(myparams$stationary_params)
  nonstationary_varnames=varname_list[!(varname_list%in%stationary_varnames)]

  for(varname in nonstationary_varnames){
    ret_augdata[[varname]]=krig_basis_to_field(
      myparams,
      pred_locs,
      varname,
      linkfuns[[varname]])
  }
  for(varname in stationary_varnames){
    ret_augdata[[varname]]=myparams$stationary_params[[varname]]
  }
  if('theta'%in%varname_list){
    ret_augdata$theta_lat=ret_augdata$theta
    ret_augdata$theta_lon=ret_augdata$theta
  }

  if(all(c('years','mu0','slope') %in% names(ret_augdata))){
    ret_augdata$mu=ret_augdata$mu0+
      (ret_augdata$years-2007)*ret_augdata$slope
  }else if('mu0' %in% names(ret_augdata)){
    ret_augdata$mu=ret_augdata$mu0
  }

  return(ret_augdata)
}

#' Parallel version of augment_data
#' @description Adds parameters to a given data frame by kriging the basis values specified in
#' "myparams" to the locations in "mydata"
#'
#' @param pred_locs Locations for kriging the parameter fields
#' @param myparams Baosis parameter values for kriging
#' @param varname_list List of parameter fields to krig
#' @param linkfuns List that contains an entry for each parameter field
#' specified in varname_list
#' @param ncores If ncores>1 parallelization is used
#' @return Returns a data frame that is a copy of the input data but with new
#' parameter values for each of the fields specified in varname_list
#' @export
augment_data_parallel=function(pred_locs,myparams,
                      varname_list=NULL,
                      linkfuns=default_linkfuns,
                      ncores=1){
  if(is.null(varname_list)){
    varname_list=names(myparams$hyperparam_list)
  }
  ret_augdata=pred_locs
  stationary_varnames=names(myparams$stationary_params)
  nonstationary_varnames=varname_list[!(varname_list%in%stationary_varnames)]

  ret_augdata[varname_list]=mclapply(varname_list,function(varname)
    krig_basis_to_field(myparams,pred_locs,varname,linkfuns[[varname]]),
    mc.cores=ncores)

  for(varname in stationary_varnames){
    ret_augdata[[varname]]=myparams$stationary_params[[varname]]
  }

  if('theta'%in%varname_list){
    ret_augdata$theta_lat=ret_augdata$theta
    ret_augdata$theta_lon=ret_augdata$theta
  }

  if(all(c('years','mu0','slope') %in% names(ret_augdata))){
    ret_augdata$mu=ret_augdata$mu0+
      (ret_augdata$years-2007)*ret_augdata$slope
  }else{
    ret_augdata$mu=ret_augdata$mu0
  }

  return(ret_augdata)
}


