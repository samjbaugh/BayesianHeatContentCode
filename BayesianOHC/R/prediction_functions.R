#'Predictions using veccmat input
#'@description Computes predictions using a cholesky precision matrix for a
#'vecchia process. Inspired by similar code in GPVecchia package; see source
#'code for citation.
#'
#'@param augdata_ordered_full Data frame containing both observations and locations
#'as well as prediction locations
#'@param pred_indices Specifies which indices of augdata_ord_full are to be predicted
#'@param veccmat The cholesky factor of the precision matrix for the Vecchia
#'process to be predicted
#'@param obs_indices The indices in augdata_ordered_full that correspond to
#'observations. By default this is the complement of the indices in pred_indices
#'@param ret_var T/F, it T calculates and returns the kriging variance
#'@param scalarvec If specified will compute predictions and variances for the
#'linear combination scalarvec*preds
#'@return If ret_var=F, returns the vector of predicted values at the specified
#'prediction locations. Otherwise returns a list with entries "pred" and
#'"var". If scalarvec is specified returns the prediction and variance for the
#'linear combination.
#'@importFrom Matrix crossprod chol solve
#'@export
pred_with_veccmat<-function(augdata_ordered_full,pred_indices,
                            veccmat,obs_indices=!pred_indices,
                      ret_var=F,scalarvec=NA){
  #this function is based on code with similar functionality in the GPvecchia
  #package by Katzfuss, Jurek, Zilber, and Gong as can be found here:
  #https://github.com/katzfuss-group/GPvecchia/ (see paper for full citation)
  all_indices=pred_indices | obs_indices
  augdata_obs=augdata_ordered_full[obs_indices,]
  zbar=(augdata_obs$z-augdata_obs$mu)/sqrt(augdata_obs$phi)
  veccmat1=veccmat[obs_indices,all_indices,drop=FALSE]
  veccmat2=veccmat[pred_indices,all_indices,drop=FALSE]
  z1=Matrix::crossprod(veccmat1,zbar)
  z2=Matrix::crossprod(Matrix::t(veccmat2),z1)
  W=Matrix::tcrossprod(veccmat2,veccmat2)
  V=Matrix::chol(W)
  z3=-Matrix::solve(Matrix::t(V),z2)
  zpred=Matrix::solve(V,z3)*
    sqrt(augdata_ordered_full$phi[pred_indices])+
    augdata_ordered_full$mu[pred_indices]
  if(ret_var){
    if(!is.na(scalarvec)){
      pred=sum(zpred*scalarvec)
      myvar=(sum(Matrix::solve(Matrix::t(V),scalarvec)^2))
    }else{
      pred=as.vector(zpred)
      tmp=Matrix::solve(Matrix::t(V),diag(sum(pred_indices)))
      myvar=Matrix::diag(Matrix::t(tmp)%*%tmp)
    }
    return(list(pred=pred,var=myvar))
  }else{
    return(as.vector(zpred))
  }
}
