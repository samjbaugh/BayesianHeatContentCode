#'Compute kriging distributions for linear combinations
#'@description Computes posterior distribution for a linear combination of
#'predicted values
#'
#'@param obsdata Data containing observations
#'@param pred_locs Locations for predictions
#'@param myparams Parameter object for kriging
#'@param grouping_list_preddf List of group information for each year
#'@param varname_list List of varnames to use
#'@param scalarvec Vector of scalars for linear combination. If missing, this
#'will be inferred by assuming pred_locs is a grid
#'@param yearlist Years to use in kriging
#'@param ncores If ncores>1 parallelization is used
#'@param ret_var T/F, if T returns the kriging variance of the linear combination
#'as well as the mean
#'@return Returns the value of pred%*%scalarvec, where pred is the vector
#'of kriging predictions located at pred_locs
compute_lincomb_predictions<-function(obsdata,pred_locs,
                                      myparams,grouping_list_preddf,
                                      varname_list=default_varname_list,
                                      scalarvec=NULL,yearlist=2007:2016,
                                   ncores=1,ret_var=T){
  if(is.null(scalarvec)){
    scalarvec=gen_deltas_from_grid(pred_locs)
  }
  obsdata$obs=T
  aug_grid=augment_data(pred_locs,myparams,varname_list)
  pred_df=obsdata[,c(names(aug_grid),'z')]
  for (year in yearlist) {
    aug_grid$years=year
    aug_grid$z=NA
    aug_grid$obs=F
    pred_df=rbind(pred_df,aug_grid)
  }

  new_veccmat_list=build_veccmat_list_grouped(pred_df,grouping_list_preddf,verb = F,ncores=ncores)
  preds_year<-function(year){
    year_df=pred_df[pred_df$years==year,]
    outobj=pred_with_veccmat(year_df,!year_df$obs,new_veccmat_list[[year-2007+1]],
                       ret_var = ret_var,scalarvec=scalarvec)
    var=ifelse(ret_var,outobj$var,NA)
    pred=ifelse(ret_var,outobj$pred,sum(outobj*scalarvec))
    return(data.frame(year=year,
                      resid=pred-sum(scalarvec*aug_grid$mu),
                      mu0=sum(scalarvec*aug_grid$mu0),
                      slope=sum(scalarvec*aug_grid$slope),
                      var=var))
  }
  predictions_ret=do.call(rbind,mclapply(2007:2016,preds_year,mc.cores = ncores))
  return(predictions_ret)
}
