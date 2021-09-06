#' Function to convert BayesNSGP samples to cv_params object
#' @param nsgp_vals Samples obtained from fitting a BayesNSGP model
#' @param knot_points Knot points used in model fit
#' @param argo_data_chord Argo data with 3D euclidean coordinates
#' @param grouping_list List of groups to use
#' @export
convert_nsgp_to_cvparams<-function(nsgp_vals,knot_points,argo_data_chord,
                                   grouping_list){
  nsgp_hypers=list(theta=list(mu=nsgp_vals$SigmaGP_mu.1.,
                              range=nsgp_vals$SigmaGP_phi.1.,
                              phi=nsgp_vals$SigmaGP_sigma,
                              nugget=1e-10),
                nugget=list(mu=nsgp_vals$tauGP_mu,
                            range=nsgp_vals$tauGP_phi,
                            phi=nsgp_vals$tauGP_sigma,
                            nugget=1e-10),
                phi=list(mu=nsgp_vals$sigmaGP_mu,
                         range=nsgp_vals$sigmaGP_phi,
                         phi=nsgp_vals$sigmaGP_sigma,
                         nugget=1e-10))

  basis_vals=list(theta=unlist(nsgp_vals[sapply(names(nsgp_vals),
                                                   function(x) substr(x,0,8)=='w1_Sigma')]),
                  nugget=unlist(nsgp_vals[sapply(names(nsgp_vals),
                                                    function(x) substr(x,0,5)=='w_tau')]),
                  phi=unlist(nsgp_vals[sapply(names(nsgp_vals),
                                                 function(x) substr(x,0,7)=='w_sigma')]))


  nsgp_basis_fields=knot_points
  nsgp_basis_fields$basis_theta=basis_vals$theta
  nsgp_basis_fields$basis_nugget=basis_vals$nugget
  nsgp_basis_fields$basis_phi=basis_vals$phi
  myparams=list(basis_fields=nsgp_basis_fields,
                hyperparam_list=nsgp_hypers)

  augdata_nsgp=argo_data_chord
  for(varname in c('phi','theta','nugget')){
    augdata_nsgp[[varname]]=krig_basis_to_field(predlocs=argo_data_chord,
                                                myparams=myparams,
                                                varname,default_linkfuns[[varname]],
                                                corrfun=chord_cor_exp_double)
  }

  betas=unlist(nsgp_vals[sapply(names(nsgp_vals),
                                function(x) substr(x,0,4)=='beta')])
  X=model.matrix(~lat_degrees*lon_degrees+
                   I(lat_degrees^2)+I(lon_degrees^2),data=augdata_nsgp)
  augdata_nsgp$mu=X%*%betas


  cv_params=build_veccmat_list_grouped(augdata_nsgp,grouping_list,
                                       corrfun=chord_cor_exp_single)

  return(list(augdata=augdata_nsgp,
              cv_params=cv_params))

}
