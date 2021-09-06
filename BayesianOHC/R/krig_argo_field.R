#' Computes kriging values and optionally variances
#'
#' @param input_data Data containing both observation locations and
#' prediction locations
#' @param pred_indices Indices of locations in the input_data object
#' on which predictions should be computed; it is implied that indicies
#' not in this vector are observation locations.
#' @param ret_var T/F, if T compute and return the kriging variane
#' @return Returns a list with entries "val" containing the kriged
#' values at the prediction locations, and "var" containing the kriging
#' variances (only if ret_var=T)
#' @export
krig_argo_field<-function(input_data,pred_indices,ret_var=F)
{
  obs_field=input_data[!pred_indices,]
  pred_field=input_data[pred_indices,]
  nobs=dim(obs_field)[1]
  npred=dim(pred_field)[1]

  c_obs_obs_mat=outer(1:nobs,1:nobs,cyl_cor_single,obs_field)
  c_obs_pred_vec=outer(1:nobs,1:npred,cyl_cor_double,obs_field,pred_field)
  myweights=solve(c_obs_obs_mat,c_obs_pred_vec)
  zbar=(obs_field$z-obs_field$mu)/sqrt(obs_field$phi)
  interp_vals=t(myweights)%*%zbar*
    sqrt(pred_field$phi)+pred_field$mu

  krig_var=NA
  if(ret_var){
    c_pred_pred=outer(1:npred,1:npred,cyl_cor_single,pred_field)
    term2=t(c_obs_pred_vec)%*%myweights
    vars=c_pred_pred-term2
  }
  return(list('pred'=interp_vals,'var'=vars))
}

#' Creates parameter field from basis values
#'
#' @param myparams List containing knot locations and basis values
#' @param predlocs Data frame containing locations for kriging
#' @param varname Name of parameter field to krig
#' @param linkfun Link function for field
#' @param corrfun Correlation function to use
#' @return Returns a vector of parameters at the locations specificed
#' by predlocs
#' @export
krig_basis_to_field<-function(myparams,predlocs,varname,linkfun,
                              corrfun=cyl_cor_double)
{
  if(varname %in% c('mu0','slope')){
    mybasis=myparams$basis_mu
  }else{
    mybasis=myparams$basis_fields
  }
  hypers=myparams$hyperparam_list
  myrange=hypers[[varname]][['range']]
  mu=hypers[[varname]][['mu']]
  phi=hypers[[varname]][['phi']]

  nobs=dim(mybasis)[1]
  npred=dim(predlocs)[1]

  predlocs$theta_lat=myrange
  predlocs$theta_lon=myrange
  predlocs$theta=myrange
  mybasis$theta_lat=myrange
  mybasis$theta_lon=myrange
  mybasis$theta=myrange

  C_pred_obs=outer(1:npred,1:nobs,corrfun,predlocs,mybasis)
  C_obs_obs=outer(1:nobs,1:nobs,corrfun,mybasis,mybasis)
  chol_obs=chol(C_obs_obs+1e-7*diag(nobs))

  varname_basis=paste('basis',varname,sep='_')
  retval=c(linkfun(sqrt(phi)*
                     C_pred_obs%*%
                     backsolve(chol_obs,mybasis[[varname_basis]])+mu))
  return(retval)
}
