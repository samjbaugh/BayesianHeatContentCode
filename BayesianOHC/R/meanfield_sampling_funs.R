#'Sample mean and trend fields
#'@description Samples the mean and trend parameters from their joint posterior
#'distribution conditioned on the other parameter fields
#'
#'@param augdata_ordered The data from which to compute the posterior.
#'@param myparams The basis values for the other parameter fields;
#'also needs to contain hyperpriors for the mean and trend fields.
#'@param yearlist List of years over which to compute the posterior
#'@param veccmat_list List of cholesky factors of precision matrices
#'for each year in yearlist.
#'@param return_posterior_mean T/F, if T returns the posterior mean, if not
#'returns a sample from the posterior distribution.
#'@param return_full_posterior T/F, if T returns the covariance matrix between
#'the mean and trend fields
#'@param ncores If ncores>1, parallelization is used
#'@param seed Default is 'rand' which uses a random seed, otherwise
#'sets a user-defined seed
#'@return If return_full_posterior=F returns a parameter object with the updated
#'basis values for the mean and slope fields. If return_full_posterior=T, returns a
#'list containing the new parameters as well as the mean and covariance
#'matrix for the posterior distribution.
#'@export
sample_mean_trend<-function(augdata_ordered,myparams,veccmat_list,
                    yearlist=c(2007:2016),return_posterior_mean=T,
                    return_full_posterior=F,ncores=1,seed='rand')
{
  basis_mu0=myparams$basis_mu
  basis_slope=myparams$basis_mu

  params_mu0=myparams$hyperparam_list$mu0
  hyperrange_mu0=params_mu0$range
  hyperphi_mu0=params_mu0$phi

  params_slope=myparams$hyperparam_list$slope
  hyperrange_slope=params_slope$range
  hyperphi_slope=params_slope$phi

  nbasis_mu0=dim(basis_mu0)[1]
  basis_mu0$theta_lat=hyperrange_mu0
  basis_mu0$theta_lon=hyperrange_mu0
  Cmu0=outer(1:nbasis_mu0,1:nbasis_mu0,cyl_cor_double,
             basis_mu0,basis_mu0)
  chol_mu0=chol(Cmu0+1e-7*diag(nbasis_mu0))

  nbasis_slope=dim(basis_slope)[1]
  basis_slope$theta_lat=hyperrange_slope
  basis_slope$theta_lon=hyperrange_slope
  Cslope=outer(1:nbasis_slope,1:nbasis_slope,cyl_cor_double,
               basis_slope,basis_slope)
  chol_slope=chol(Cslope+1e-7*diag(nbasis_slope))

  compute_terms_by_year<-function(year){
    predlocs_ordered=augdata_ordered[augdata_ordered$years==year,]
    npred=dim(predlocs_ordered)[1]
    myphi_ordered=predlocs_ordered$phi
    veccmat=veccmat_list[[year-2007+1]]

    ####Compute A
    predlocs_ordered$theta_lat=hyperrange_mu0
    predlocs_ordered$theta_lon=hyperrange_mu0
    C_pred_basis_mu0=outer(1:npred,1:nbasis_mu0,cyl_cor_double,
                           predlocs_ordered,basis_mu0)

    A_ordered=sqrt(hyperphi_mu0) * C_pred_basis_mu0 %*% solve(chol_mu0)

    ####Compute B
    predlocs_ordered$theta_lat=hyperrange_slope
    predlocs_ordered$theta_lon=hyperrange_slope
    C_pred_basis=outer(1:npred,1:nbasis_slope,cyl_cor_double,
                       predlocs_ordered,basis_slope)

    B_ordered=sqrt(hyperphi_slope) * C_pred_basis %*% solve(chol_slope)


    A_aug=cbind(rep(1,npred),(year-2007)*B_ordered,
                A_ordered)
    Astand=apply(A_aug,2,function(x) x/sqrt(myphi_ordered))

    ####Done Compute A
    zcnorm_ord=(predlocs_ordered$z)/sqrt(myphi_ordered)
    zvec=Matrix::crossprod(veccmat,zcnorm_ord)
    NumA=t(Astand)%*%Matrix::crossprod(Matrix::t(veccmat),zvec)

    temp=Matrix::crossprod(veccmat,Astand)
    DenomA=Matrix::t(temp)%*%temp

    return(list('NumA'=NumA,'DenomA'=DenomA))
  }

  A_out=mclapply(yearlist,compute_terms_by_year,mc.cores = ncores)
  outvec=solve(Reduce('+',lapply(A_out,function(x) x$DenomA)),
               Reduce('+',lapply(A_out,function(x) x$NumA)))

  slope_inds=2:(nbasis_slope+1)
  mu0_inds=(nbasis_slope+2):(nbasis_slope+nbasis_mu0+1)

  if(!return_posterior_mean){
    if(seed!='rand'){
      set.seed(seed)
    }
    bmufull_var=solve(Reduce('+',lapply(A_out,function(x) x$DenomA)))
    bmufull_var_sym=(bmufull_var+t(bmufull_var))/2
    bmufull_sample=as.vector(outvec +
                               t(chol(bmufull_var_sym))%*%rnorm(1+nbasis_mu0+nbasis_slope))

    newc=bmufull_sample[1]
    bslope_ret=bmufull_sample[slope_inds]
    bmu_ret=bmufull_sample[mu0_inds]
  }else{
    newc=outvec[1]
    bslope_ret=outvec[slope_inds]
    bmu_ret=outvec[mu0_inds]
  }

  # ordered_res=A_ordered %*% bmu_mean+c
  retparams=myparams
  retparams$basis_mu$basis_mu0=bmu_ret
  retparams$basis_mu$basis_slope=bslope_ret
  retparams$hyperparam_list$mu0$mu=newc
  retparams$hyperparam_list$slope$mu=0
  if(return_full_posterior){
    bmufull_var=solve(Reduce('+',lapply(A_out,function(x) x$DenomA)))
    bmufull_mean=outvec
    return(list('retparams'=retparams,
                'mean'=bmufull_mean,
                'var'=bmufull_var))
  }else{
    return(retparams)
  }
}

#'Computes mean/trend integrated posteriors
#'@description Computes the mean and variance for the posterior distribution of
#'the globally integrated mean and trend fields as well as their
#'correlation
#'
#'@param sampleout The output of sample_meanfield
#'@param myparams The current iteration's parameter object
#'@param pred_locs Grid for computing integrated field
#'@return Returns a list containing mu0 (the posterior mean of the
#'integrated mean field), mu0_var (the posterior variance of the
#'integrated mean field), slope and slope_var analogously, and
#'corr_term which gives the correlation between the integrated
#'mean and slopes
#'@export
get_mu_integrated_posterior<-function(sampleout,myparams,pred_locs){
  myrange=myparams$hyperparam_list$slope$range

  mybasis=myparams$basis_mu
  mybasis$theta_lon=myrange
  mybasis$theta_lat=myrange
  mybasis$nugget=1e-7

  scalarvec=gen_deltas_from_grid(pred_locs)

  mypreds=pred_locs
  mypreds$theta_lon=myrange
  mypreds$theta_lat=myrange
  nbasis=dim(mybasis)[1]
  nnew=dim(mypreds)[1]

  C_new_old=outer(1:nnew,1:nbasis,cyl_cor_double,mypreds,mybasis)
  C_old_old=outer(1:nbasis,1:nbasis,cyl_cor_single,mybasis)
  chol_old=chol(C_old_old)

  weight_matrix_c=rep(1,nnew)
  weight_matrix_mu=sqrt(myparams$hyperparam_list$mu0$phi)*C_new_old%*%solve(chol_old)
  weight_matrix_slope=sqrt(myparams$hyperparam_list$slope$phi)*C_new_old%*%solve(chol_old)
  weight_matrix_total=cbind(weight_matrix_c,weight_matrix_slope,weight_matrix_mu)

  mu0_inds=c(1,(nbasis+2):(nbasis*2+1))
  slope_inds=2:(nbasis+1)
  tempfun<-function(aa,bb){aa[-bb]=0;return(aa)}
  onemat=cbind(tempfun(rep(1,nbasis*2+1),mu0_inds),
               tempfun(rep(1,nbasis*2+1),slope_inds))
  all_means=(scalarvec%*%weight_matrix_total%*%diag(c(as.matrix(sampleout$mean)))%*%onemat)
  covmat=t(onemat)%*%diag(c(scalarvec%*%weight_matrix_total))%*%sampleout$var%*%
    diag(c(scalarvec%*%weight_matrix_total))%*%onemat
  covmat_sym=(covmat+Matrix::t(covmat))/2
  mycor=covmat_sym[1,2]/(sqrt(covmat_sym[1,1]*covmat_sym[2,2]))
  retval=data.frame('mu0'=all_means[1],'mu0_var'=sqrt(covmat[1,1]),
                    'slope'=all_means[2],'slope_var'=sqrt(covmat[2,2]),
                    'corr_term'=mycor)
  retval$zscore=retval$slope/retval$slope_var
  return(retval)
}
