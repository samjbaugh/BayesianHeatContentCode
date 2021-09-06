#' Re-samples slopes from marginal posterior
#' @description Re-samples slopes from the marginal posterior conditional on
#' the values of the other parameters.
#'
#' @param mu_sampleout List containing parameters and posterior covariances
#' for the mean and trend fields. Should be the output of "sample_mean_trend"
#' run with return_full_posterior=T.
#' @param nresamp Number of re-sampled fields to return.
#' @param krig_grid Grid for computing integrated mean/trend values.
#' @param scalarvec Values for numerical integrated; if missing will infer
#' from grid.
#' @return Returns a list with "trend_samples", a nresample(x)ngrid matrix
#' where each row contains the trend values of a resampled field at the
#' corresponding grid-points, "intslope_samples" which is a vector of length
#' nresamp giving the integrated values of "trend_samples", "intmean" which
#' gives the posterior mean of the integrated trend, and "intvar" which
#' gives the posterior variance of the integrated trend.
#' @export
resample_slopes<-function(mu_sampleout,nresamp,krig_grid,
                          scalarvec=NULL){
  ###get stuff section
  if(is.null(scalarvec)){
    scalarvec=gen_deltas_from_grid(krig_grid)
  }
  myparams=mu_sampleout$retparams
  myrange=myparams$hyperparam_list$slope$range

  mybasis=myparams$basis_mu
  mybasis$theta_lon=myrange
  mybasis$theta_lat=myrange
  mybasis$nugget=1e-7

  mypreds=krig_grid
  mypreds$theta_lon=myrange
  mypreds$theta_lat=myrange
  nold=dim(mybasis)[1]
  nnew=dim(mypreds)[1]

  C_new_old=outer(1:nnew,1:nold,cyl_cor_double,mypreds,mybasis)
  C_old_old=outer(1:nold,1:nold,cyl_cor_single,mybasis)
  chol_old=chol(C_old_old)

  weight_matrix_c=rep(1,nnew)
  weightmatbase=C_new_old%*%solve(chol_old)
  weight_matrix_mu=sqrt(myparams$hyperparam_list$mu0$phi)*weightmatbase
  weight_matrix_slope=sqrt(myparams$hyperparam_list$slope$phi)*weightmatbase
  weight_matrix_total=cbind(weight_matrix_c,weight_matrix_slope,weight_matrix_mu)
  mu0_inds=c(1,260:517)
  slope_inds=2:259

  onemat=cbind(rep(1,517),rep(1,517))
  onemat[slope_inds,1]=0
  onemat[mu0_inds,2]=0

  slope_meanfield=(weight_matrix_total%*%(diag(c(as.matrix(mu_sampleout$mean)))%*%onemat))[,2]
  intslope_mean=sum(scalarvec*slope_meanfield)
  varmat=t(onemat)%*%diag(c(scalarvec%*%weight_matrix_total))%*%mu_sampleout$var%*%
    diag(c(scalarvec%*%weight_matrix_total))%*%onemat
  intslope_var=sqrt(varmat[2,2])

  #only for numerical reasons:
  basis_covmat=(t(as.matrix(mu_sampleout$var))+as.matrix(mu_sampleout$var))/2
  chol_covmat=chol(basis_covmat)

  randmat=matrix(rnorm(517*nresamp),517,nresamp)
  newsample_basis=mu_sampleout$mean+t(chol_covmat)%*%randmat
  trend_samples=t(as.matrix(weight_matrix_total%*%(newsample_basis*onemat[,2])))
  intslope_samples=apply(trend_samples,1,function(x) sum(scalarvec*x))

  retval=list('trend_samples'=trend_samples,'intslope_samples'=intslope_samples,
              'intmean'=intslope_mean,'intvar'=intslope_var)
  return(retval)
}
