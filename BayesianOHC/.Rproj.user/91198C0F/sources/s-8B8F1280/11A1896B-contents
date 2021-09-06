#' Cholesky likelihood function
#'
#' @param z Simulated observed values
#' @param ranges Vector of range values
#' @param locs Locations
#' @param corfun Function for computing correlations
#' @export
chol_lik<-function(z,ranges,locs,corfun){
  nlocs=length(locs)
  nuggetparam=1e-3
  retval=tryCatch({
    nz=length(z)
    covmat=outer(1:nlocs,1:nlocs,function(ii,jj)
      corfun(locs[ii],locs[jj],ranges[ii],ranges[jj]))+
      nuggetparam*diag(nlocs)
    mychol=chol(covmat)

    quadterm=-sum(apply(forwardsolve(t(mychol),z)^2,2,sum))
    diagterm=-2*sum(log(diag(mychol)))#-npixel*log(scale)
    nz=dim(z)[1]
    piterm=-nz*log(2*pi)
    retval=sum(quadterm+diagterm+piterm)/2
    retval
  },error=function(e){print(e); print(max(ranges)); -Inf})
  return(retval)
}

#' Simulates data using true parameters, and finds MLEs
#'
#' @description MLEs are found using both cylindrical_correlation_exact anad
#'  cylindrical_correlation_exact
#' @param true_ranges The true ranges to use for the simulation
#' @param rangefun Function for computing non-stationary range parameters
#' @param locs Locations on which to generate data
#' @export
simulate_mles<-function(true_ranges,rangefun,locs){
  nlocs=length(locs)
  nuggetparam=1e-3
  # locs=runif(100,-pi,pi)
  covmat=outer(1:nlocs,1:nlocs,function(ii,jj)
    cylindrical_correlation_exact(locs[ii],locs[jj],true_ranges[ii],true_ranges[jj]))+
    nuggetparam*diag(nlocs)
  simdata=t(chol(covmat))%*%rnorm(nlocs)

  optfun_conv=function(param) -chol_lik(simdata,rangefun(param),
                                        corfun=cylindrical_correlation_exact,
                                        locs=locs)
  optfun_approx=function(param) -chol_lik(simdata,rangefun(param),
                                          corfun=cylindrical_correlation_gaussian,
                                          locs=locs)

  optout_conv=optimize(optfun_conv,lower=0,upper=2*max(true_ranges))
  optout_approx=optimize(optfun_approx,lower=0,upper=2*max(true_ranges))

  x=list(conv=optout_conv$minimum,
         approx=optout_approx$minimum)
  return(x)
}

#' Runs simulation comparing exact cylindrical convolutions with
#' Gaussian approximation
#'
#' @param true_range_seq Sequence of true range parameters
#' @param rangefun Function for computing non-stationary range parameters
#' @param nrep Number of repetitions to simulate data
#' @param locs Locations on which to generate data
#' @export
cylind_approx_simulation<-function(true_range_seq,rangefun,nrep,locs){
  nlocs=length(locs)
  nranges=length(true_range_seq)
  mle_conv=matrix(NA,nranges,nrep)
  mle_approx=matrix(NA,nranges,nrep)
  for(aii in 1:nranges){
    true_range=true_range_seq[aii]
    true_ranges=rangefun(true_range)

    print(paste('aii:',aii,' true_range:',true_range,sep=''))

    opt_wrapper<-function(x) simulate_mles(true_ranges,rangefun,locs)

    aa=mclapply(1:nrep,function(ii) opt_wrapper(),mc.cores=8)

    mle_approx[aii,]=sapply(aa,function(x) x$approx)
    mle_conv[aii,]=sapply(aa,function(x) x$conv)
  }
  simdf=data.frame(true_range_seq=true_range_seq,
               conv=apply(mle_conv,1,mean),
               approx=apply(mle_approx,1,mean),
               true_effrange=convert_theta_lon_to_effective_range_deg(true_range_seq))
  simdf$conv_effrange=convert_theta_lon_to_effective_range_deg(simdf$conv)
  simdf$approx_effrange=convert_theta_lon_to_effective_range_deg(simdf$approx)
  return(simdf)
}
