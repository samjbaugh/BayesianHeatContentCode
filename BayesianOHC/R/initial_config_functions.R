#'Find MLEs for stationary process
#' @description Computes maximum likelihood estimates for an isotropic cylindrical model
#'
#' @param mydata dataframe with column z containing osbervations
#' @param yearlist list of years to use, default is 2007:2016
#' @param verb If T use trace=6 in optim
#' @param prior_list The list of coefficients for normal priors
#' on the log variables (if NULL or missing no priors are used)
#' @param est_nugget T/F if T nugget is estimated
#' @param est_mu T/F, if T constant mean is estimated
#' @param est_slope T/F, if T linear slope is estimated
#' @param iso T/F, if T an isotropic model will be used
#' @return Returns a list containing the maximum likelihood estimates for
#' a stationary cylindrical model on the given data. Output contains
#' "theta_lat", "theta_lon", "nugget", "mu0", "phi", "slope", "likelihood",
#' "data_mean" and "data_var"
#' @import stats
#' @export
get_stationary_MLEs <- function(mydata,yearlist=2007:2016,verb=F,prior_list=NULL,
                            est_nugget=F,est_mu=T,est_slope=T,iso=F)
{
  myn=dim(mydata)[1]
  #non-isotropic case
  start_value=c(-5,5)
  if(est_nugget)
  {
    start_value=c(start_value,-5)
  }else{
    fixed_nugget=1e-2
  }

  optfun=function(myparams,return_profiled=F) {
    inputparams=c(myparams[1],ifelse(iso,myparams[1],myparams[2]),
                  ifelse(est_nugget,myparams[3],log(fixed_nugget)))
    retval = tryCatch({
      gp_profiled_stationary_likelihood(mydata,inputparams,
                                        yearlist=yearlist,
                                        prior_list=prior_list,return_profiled=return_profiled,
                                        est_mu=est_mu,est_slope=est_slope,
                                        iso=iso)
    }, warning = function(w) {
      return(-9e16)
    }, error = function(e) {
      return(-9e16)
    })
    return(retval)
  }
  trace=ifelse(verb,6,0)

  if(iso & !est_nugget){
    optobj=optimize(function(x) -optfun(x),lower=-20,upper=10) #
    gpout=optfun(optobj$minimum,return_profiled=T)
    retval=list("range"=exp(optobj$minimum),"nugget"=fixed_nugget,
              "mu0"=gpout$mu0,"phi"=gpout$phi,slope=gpout$slope,
              "likelihood"=-optobj$objective,"window_size"=length(mydata$z),
              "data_mean"=mean(mydata$z),"data_var"=var(mydata$z))
  }else{
    optobj=optim(par=start_value,fn=optfun,method="BFGS",
                 control=list(trace=trace,fnscale=-1)) #,control=list(maximize=T,info=T))
    gpout=optfun(optobj$par,return_profiled=T)
    retval=list("theta_lat"=exp(optobj$par)[1],"theta_lon"=exp(optobj$par)[2],
                "nugget"=ifelse(est_nugget,exp(optobj$par)[3],fixed_nugget),
                "mu0"=gpout$mu0,"phi"=gpout$phi,slope=gpout$slope,
                "likelihood"=-optobj$value,"window_size"=length(mydata$z),
                "data_mean"=mean(mydata$z),"data_var"=var(mydata$z))
  }

  # "liks"=gpout$liks,"priors"=gpout$priors)
  return(retval)
}

#'Computes profiled log-likelihood for stationary Gaussian process
#' @description Compute Gaussian process log-likelihood with options to profile
#' over the mean and variance parameters
#'
#' @param mydata dataframe with column z containing osbervations
#' @param logparams vector of log parameters (theta_lat,theta_lon,nugget)
#' @param yearlist list of years to use, default is 2007:2016
#' @param return_profiled T/F, if T profiled parameters are returned
#' @param prior_list The list of coefficients for normal priors
#' on the log variables (if NULL or missing no priors are used)
#' @param est_mu T/F, if T constant mean is estimated
#' @param est_slope T/F, if T linear slope is estimated
#' @param iso T/F, if T an isotropic model will be used
#' @return If return_profiled=F, returns the log-likelihood. Otherwise,
#' returns a list containing the log-likelihood and the profiled mu0, slope,
#' and phi (variance) values
#' @export
gp_profiled_stationary_likelihood<-function(mydata,logparams,
                               yearlist=2007:2016,
                               return_profiled=F,
                               prior_list=NULL,est_mu=F,
                               est_slope=F,iso=F)
{
  use_prior=!is.null(prior_list)

  theta_lat=exp(logparams[1])
  theta_lon=exp(logparams[2])
  nugget=exp(logparams[3])
  if(iso){
    theta_lon=theta_lat
  }

  mydata$theta_lat=theta_lat
  mydata$theta_lon=theta_lon
  mydata$nugget=nugget

  if(use_prior)
  {
    phi_prior_term=0 #-(log(phi)-prior_list$mu_phi)^2/(2*prior_list$sigma_phi^2)
    nugget_prior_term=-(log(nugget)-prior_list$mu_nugget)^2/(2*prior_list$sigma_nugget^2) #ifelse(!prior,0,-(log(nugget)-mu_nugget)^2/(2*sigma_nugget^2)) #-(alpha+1)*log(nugget)-beta/nugget
    lat_prior_term=-(log(theta_lat)-prior_list$mu_lat)^2/(2*prior_list$sigma_lat^2) #-(alpha+1)*log(nugget)-beta/nugget
    lon_prior_term=-(log(theta_lon)-prior_list$mu_lon)^2/(2*prior_list$sigma_lon^2) #-(alpha+1)*log(nugget)-beta/nugget
    if(iso){
      lon_prior_term=lat_prior_term
    }
    priors=c(lat_prior_term,lon_prior_term,nugget_prior_term)
  }else{
    priors=c(0)
  }
  #
  liks=rep(NA,length(yearlist))
  covmat_chol_years=list()

  get_mu_term<-function(year_i){
    yeardata=mydata[mydata$years==yearlist[year_i],]
    ndata=dim(yeardata)[1]
    if(ndata==0){
      return(list(covmat_chol=0,mu_numer=0,mu_denom=0))
    }
    covmat_full=outer(1:ndata,1:ndata,cyl_cor_single,yeardata)
    covmat_chol=chol(covmat_full)

    covmat_chol_years[[year_i]]=covmat_chol
    if(est_slope){
      X=cbind(rep(1,ndata),(yearlist[year_i]-2007)*rep(1,ndata))
    }else{
      X=cbind(rep(1,ndata))
    }
    mu_numer=t(X)%*%(chol_solve(covmat_chol,yeardata$z))
    temp=forwardsolve(t(covmat_chol),X)
    mu_denom=t(temp)%*%temp
    return(list(covmat_chol=covmat_chol,mu_numer=mu_numer,mu_denom=mu_denom))
  }
  outterms=lapply(1:length(yearlist),get_mu_term)
  if(est_mu & est_slope){
    outvec=solve(Reduce('+',lapply(outterms,function(x) x$mu_denom)),
                 Reduce('+',lapply(outterms,function(x) x$mu_numer)))
    mu0=outvec[1]
    slope=outvec[2]
  }else{
    if(est_mu & !est_slope){
      outvec=solve(Reduce('+',lapply(outterms,function(x) x$mu_denom)),
                   Reduce('+',lapply(outterms,function(x) x$mu_numer)))
      mu0=outvec[1]
      slope=0
    }else{
      mu0=0 #outvec[1]
      slope=0
    }
  }

  get_phiterm<-function(year_i){
    year=yearlist[year_i]
    zbar=mydata[mydata$years==yearlist[year_i],]$z-mu0-slope*(year-2007)
    covmat_chol=outterms[[year_i]]$covmat_chol
    if(length(zbar)==0){
      return(0)
    }
    return(sum(forwardsolve(t(covmat_chol),zbar)^2))
  }

  phi=sum(sapply(1:length(yearlist),get_phiterm))/sum(sapply(yearlist,function(x) sum(mydata$years==x)))

  get_likterm=function(year_i){
    year=yearlist[year_i]
    zbar=(mydata[mydata$years==yearlist[year_i],]$z-mu0-slope*(year-2007))
    if(length(zbar)==0){
      return(0)
    }
    covmat_chol=outterms[[year_i]]$covmat_chol * sqrt(phi)
    deranged_z=forwardsolve(t(covmat_chol),zbar)
    term1=2*sum(log(diag(covmat_chol)))
    term2=sum(deranged_z^2)
    return(-(term1+term2)/2)
  }
  liks=sapply(1:length(yearlist),get_likterm)

  if(return_profiled)
  {
    return(list('likelihood'=sum(liks)+sum(priors),
                "mu0"=mu0,"phi"=phi,"slope"=slope,"priors"=priors,'liks'=liks))
  }
  else
  {
    return(sum(liks)+sum(priors))
  }
}

