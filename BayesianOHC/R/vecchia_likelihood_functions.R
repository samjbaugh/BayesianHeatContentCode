#'Vecchia likelihood from veccmat object
#' @description Computes the vecchia likelihood from the cholesky of the precision
#' matrix of a vecchia process. Based off of code with similar functionality in
#' the GPVecchia package, see source code for citation details.
#'
#' @param likinput_ordered The ordered data on which to evaluate the likelihood,
#' where observational values are stored in the z column
#' @param veccmat The cholesky of the precision matrix corresponding to a
#' vecchia process on the data
#' @return Returns the log-likelihood of the observations
#' @export
vecc_lik_from_veccmat<-function(likinput_ordered,veccmat) {
  #this function is inspired by code with similar functionality in the GPvecchia
  #package by Katzfuss, Jurek, Zilber, and Gong as can be found here:
  #https://github.com/katzfuss-group/GPvecchia/ (see paper for full citation)
  zord=(likinput_ordered$z-likinput_ordered$mu)/sqrt(likinput_ordered$phi)
  const=dim(likinput_ordered)[1]*log(2*pi)
  z=Matrix::crossprod(veccmat,zord)
  quadterm=sum(z^2)
  logdetterm=-2*sum(log(Matrix::diag(veccmat)))
  phiterm=sum(log(likinput_ordered$phi))
  neg2loglik=logdetterm+quadterm+const+phiterm
  loglik=-neg2loglik/2
  return(loglik)
}

#' Evaluates the vecchia likelihood over multiple years
#'
#' @param likinput_ordered The ordered data on which to evaluate the likelihood,
#' where observational values are stored in the z column
#' @param veccmat_list List of cholesky of the precision matrices corresponding
#' to the years specified in yearlist
#' @param yearlist The years over which to calculate the likelihoods
#' @return Returns a list of log-likelihoods the same length as
#' yearlist
#' @export
vecc_lik_over_years=function(likinput_ordered,
                                   veccmat_list,yearlist=2007:2016){
  lik_list=lapply(yearlist,function(myyear){
    yeardata=likinput_ordered[likinput_ordered$years==myyear,]
    vecc_lik_from_veccmat(yeardata,veccmat_list[[myyear-min(yearlist)+1]])
  })
  return(unlist(lik_list))
}
