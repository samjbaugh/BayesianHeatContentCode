#' Computes cross-validation scores and standard errors
#' @description Computes cross-validation scores and their associated prediction
#' errors using either lofo (leave-one-float-out) or lowo (leave-one-window-out)
#' cross validation methods
#'
#' @param augdata_ordered The input data augmented with the kernel parameter
#' values at each location (can be created with function "augment_data")
#' @param cv_function Default is "cv_scores", an availabel option is
#'  "levitus_preds".
#' @param cv_params List of parameters to the cross-validation function;
#' should be the same length as yearlist. If cvtype=standard this is a list of
#' veccmat objects, if cvtype=levitus this should be a distance matrix between
#' each of the observations in each year in kilometers, and if cvtype=gpgp this
#' should be a fitted GpGp model
#' @param yearlist The list of years over which to compute the validation scores
#' @param cv_type lofo (leave-one-float-out) vs lowo (leave-one-window-out)
#' @param nsamp If nsamp=-1 computes scores for all observation locations; if
#' not validation scores are computed for a selection of nsamp indices for
#' each year
#' @param cv_options Options for cv; if lowo is chosen you can input the width
#' of the window to withold
#' @param ncores If ncores>1, how many cores should be used for parallelization
#' @return Returns a dataframe the same size of "augdata_ordered" with two
#' additional columns, validval and validsd, giving the mean and standard
#' deviation of the leave-one-float-out validation score at that point
#' @export
compute_cv_df<-function(augdata_ordered,cv_params,cv_function=cv_scores,
                        cv_type='lofo',yearlist=2007:2016,nsamp=-1,
                        ncores=1,cv_options=NULL){
  validdf_list=list()
  for(myyear in yearlist){
    print(myyear)
    yeardata_ordered=augdata_ordered[augdata_ordered$years==myyear,]

    if(cv_type=='lofo'){
      float_ids=unique(yeardata_ordered$float_id)
      if(nsamp>0){
        set.seed(123)
        float_ids=sample(float_ids,nsamp)
      }
      cv_indices=lapply(float_ids,function(fid) yeardata_ordered$float_id==fid)
    }
    if(cv_type=='lowo'){
      nyear=dim(yeardata_ordered)[1]
      indexlist=1:nyear
      if(nsamp>0){
        indexlist=sample(indexlist,nsamp)
      }
      cv_indices=lapply(indexlist,function(obsid)
        c(window_dist(yeardata_ordered,yeardata_ordered[obsid,])<
             deg2rad(cv_options$lowo_width/2)))

    }
    validlist=mclapply(cv_indices,cv_function,yeardata_ordered,
                       cv_params[[myyear-2007+1]],mc.cores=ncores)

    validdf=data.frame(do.call(rbind,validlist))
    validdf_list[[myyear-2007+1]]=validdf
  }
  retval=do.call(rbind,validdf_list)
  return(retval)
}

compute_cv_df_wrapper<-function(augdata_ordered,
                                cv_type='lofo',cv_options=NULL,
                                yearlist=2007:2016,nsamp=-1,
                                ncores=1){
  grouping_list=compute_grouping_list(augdata_ordered,m=50,yearlist=yearlist)
  cv_params=build_veccmat_list_grouped(augdata_ordered,
                                       grouping_list,ncores=ncores,
                                       yearlist=yearlist)
  retval=compute_cv_df(augdata_ordered,cv_function=cv_scores,
                       cv_params=cv_params,cv_type=cv_type,
                       yearlist=yearlist,nsamp=nsamp,
                       cv_options=cv_options,
                       ncores=ncores)
  return(retval)
}

#'Convert GpGp model fit object to veccmat object
#'
#'@param gpgp_fit_list List of model fits GpGp::fit_model
#'@param ordered_data Observation data
#'@param yearlist List of years to use, should be same length as gpgp_fit_list
#'@param corrfun Correlation function to use
#'@param yearlist Lits of years to use
#'@return Returns a list of "veccmat" objects
#'@export
convert_gpgp_fit_to_cvparams<-function(gpgp_fit_list,ordered_data,corrfun,
                                       yearlist=2007:2016){
  cv_params=list()
  for(myyear in yearlist){
    print(myyear)
    gpgp_model_fit=gpgp_fit_list[[myyear-2007+1]]
    yeardata=ordered_data[ordered_data$years==myyear,]

    mygroupobj=compute_grouped_conditioning_sets(gpgp_model_fit$NNarray)
    myvecc_chordal=create_veccmat_grouped(yeardata,
                                          mygroupobj,
                                          corrfun=corrfun,
                                          covparms=gpgp_model_fit$covparms)
    cv_params[[myyear-2007+1]]=myvecc_chordal
  }
  return(cv_params)
}

#' Compute cross-validation scores with veccmat object
#'
#' @param cv_indices Indices to predict in cross-validation
#' @param obspred_df Dataframe containing both observation and prediction
#' locations
#' @param veccmat Veccmat object corresponding to data in obspred_df
#' @export
cv_scores<-function(cv_indices,obspred_df,veccmat)
{
  validout=pred_with_veccmat(obspred_df,cv_indices,
                             veccmat,ret_var = T)
  validdf=obspred_df[cv_indices,]
  validdf$validval=validout$pred
  validdf$validsd=sqrt(validout$var)
  return(validdf)
}

#' Compute predictions using method described by Levitus et al. (2012),
#' see supplementary materials for details
#'
#' @param pred_indices Indices to predict
#' @param obspred_df Dataframe containing observation and prediction locations
#' @param distmat_km Optional, can make it faster if doing cross validation
#' than having to re-compute distances
#' @param distfun If distmat not supplied, function for computing distances;
#' should be geodist for great circle distances
#' @param verb Should progress be displayed?
#' @return Dataframe with "validval" (predictions) and "validsd"
#' (standard errors)
#' @export
levitus_preds<-function(pred_indices,obspred_df,distmat_km=NULL,
                        distfun=NULL,verb=F){
  RR=666
  E=4
  validdf=obspred_df[pred_indices,,drop=F]
  obsdf=obspred_df[!pred_indices,,drop=F]
  calc_dists=is.null(distmat_km)
  if(!calc_dists){
    mydists=distmat_km[!pred_indices,pred_indices,drop=F]
  }
  for(predi in 1:sum(pred_indices)){
    if(verb){
      print(predi)
    }
    if(calc_dists){
      mydists=distfun(obspred_df[!pred_indices,c('lon_degrees','lat_degrees')],
                      obspred_df[which(pred_indices)[predi],c('lon_degrees','lat_degrees')],
                      measure='haversine')/1e3
      r=mydists
    }else{
      r=mydists[,predi,drop=F]
    }
    in_region=r<RR
    weights=exp(-E*r[in_region]^2/RR^2)
    validdf$validval[predi]=sum(weights*obsdf$z[in_region])/
      sum(weights)+validdf$mu0[predi]
    sigma0=sd(obsdf$vhc_obs[in_region])
    validdf$validsd[predi]=(sigma0/sum(weights))*
      sqrt(sum(weights%*%t(weights)))
    #assume zero anomaly if no information:
    if(sum(in_region)==0){validdf$validval[predi]=validdf$mu0[predi]}
    #if no sd information use first pass value
    if(is.na(sigma0)){validdf$validsd[predi]=validdf$sd[predi]}
  }
  return(validdf)
}

#' Add means from GpGp model fit back to data
#' @param yeardata Data from a particular year
#' @param model_fit GpGp model fit containing design matrix and parameters
#' @export
add_means<-function(yeardata,model_fit){
  covparms=model_fit$covparms
  beta=model_fit$betahat
  X=model_fit$X

  yeardata$mu=c(X%*%beta)
  yeardata$phi=covparms[1]
  return(yeardata)
}
