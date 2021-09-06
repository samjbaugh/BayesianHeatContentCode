#' Restricts the dataframe given to a specified region.
#'
#' @param df Dataframe containing locations in columns "lon_degrees" and
#' "lat_degrees".
#' @param region_name Options are North Atlantic (AtlN), Pacific Northeast (PacNE),
#' Pacific Northwest (PacNW), Tropical Atlantic (AtlTrop), Western Tropical
#' Pacific (PacTropW), Eastern Tropical Pacific (PacTropE), Indian (Ind),
#' Southern Pacific (PacSo), Southern Indian (IndSo), Southern Atlantic
#' (AtlSo), and globe. See manuscript or source code for region
#' definitions
#' @return Returns a subset of df with locations constrained to the
#' specified region.
#' @export
in_region<-function(df,region_name) {
  if(region_name=='globe'){
    return(df)
  }
  region_limits=list('AtlN'=list('latlims'=c(20,65),'lonlims'=c(-80,0)),
                     'PacNE'=list('latlims'=c(20,65),'lonlims'=c(-180,-80)),
                     'PacNW'=list('latlims'=c(20,65),'lonlims'=c(120,180)),
                     'AtlTrop'=list('latlims'=c(-30,20),'lonlims'=c(-70,30)),
                     'PacTropW'=list('latlims'=c(-30,20),'lonlims'=c(-180,-70)),
                     'PacTropE'=list('latlims'=c(-30,20),'lonlims'=c(125,180)),
                     'Ind'=list('latlims'=c(-30,20),'lonlims'=c(30,120)),
                     'PacSo'=list('latlims'=c(-80,-30),'lonlims'=c(-180,-60)),
                     'IndSo'=list('latlims'=c(-80,-30),'lonlims'=c(30,180)),
                     'AtlSo'=list('latlims'=c(-80,-30),'lonlims'=c(-60,30)),
                     'globe'=list('latlims'=c(-90,90),'lonlims'=c(-180,180)))

  lonlims=region_limits[[region_name]][['lonlims']]
  latlims=region_limits[[region_name]][['latlims']]
  lonlims2=region_limits[[region_name]][['lonlims2']]
  latlims2=region_limits[[region_name]][['latlims2']]
  return(df[((df$lon_degrees>=lonlims[1] & df$lon_degrees<=lonlims[2])) &
                  (df$lat_degrees>=latlims[1] & df$lat_degrees<=latlims[2]),])
}

#' For basis parameter, restricts the knot locations to a specific region and
#' re-krigs the parameter values.
#'
#' @param inputparams Basis function parameters
#' @param region_name Region name, see documentation for "in_region" for
#' available options
#' @param varname_list List of names of parameters to convert
#' @param invlinkfuns Inverse link functions for parameter fields
#' @return Returns a list of parameters with the same hyperparameters
#' as inputparameters
#' @export
get_region_params<-function(inputparams,region_name,
                            varname_list=default_varname_list,
                            invlinkfuns=default_invlinkfuns)
{
  if(region_name=='globe'){
    return(inputparams)
  }
  regionparams=inputparams
  regionparams$basis_fields=in_region(regionparams$basis_fields,region_name)
  regionparams$basis_mu=in_region(regionparams$basis_mu,region_name)

  region_gridparams=augment_data(regionparams$basis_fields,inputparams,varname_list)
  newparams=regionparams
  for(varname in varname_list){
    myrange=inputparams$hyperparam_list[[varname]]$range
    region_gridparams$theta_lat=myrange
    region_gridparams$theta_lon=myrange

    if(varname %in% c('mu0','slope')){
      mybasis=inputparams$basis_mu
    }else{
      mybasis=inputparams$basis_fields
    }
    mybasis$theta_lat=myrange
    mybasis$theta_lon=myrange

    npred=dim(region_gridparams)[1]
    nobs=dim(mybasis)[1]

    C_pred_obs=outer(1:npred,1:nobs,cyl_cor_double,region_gridparams,mybasis)
    C_obs_obs=outer(1:nobs,1:nobs,cyl_cor_double,mybasis,mybasis)
    chol_obs=chol(C_obs_obs+1e-10*diag(nobs))

    varname_basis=paste('basis',varname,sep='_')
    retval=C_pred_obs%*%backsolve(chol_obs,mybasis[[varname_basis]])

    C_pred_pred=outer(1:npred,1:npred,cyl_cor_double,
                      region_gridparams,region_gridparams)

    invlink=invlinkfuns[[varname]]
    chol_pred=chol(C_pred_pred+1e-10*diag(npred))
    newbasis=solve(t(chol_pred),retval)

    if(varname %in% c('mu0','slope')){
      newparams$basis_mu[[paste('basis',varname,sep='_')]]=newbasis
    }else{
      newparams$basis_fields[[paste('basis',varname,sep='_')]]=newbasis
    }
  }
  return(newparams)
}
