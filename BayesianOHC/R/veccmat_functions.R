#' Computes list of cholesky factors of the Vecchia precision matrix
#'
#' @param augdata_ordered Ordered data over which to compute the factors. Data
#' should be 'augmented' meaning that each row has the values of the kernel
#' parametmers for that locations; this can be added to a dataframe by calling
#' "augment_data" with a parameter object
#' @param grouping_list The list of grouping information lists for
#' the process
#' @param yearlist The list of years for creating the factors
#' @param corrfun Correlation function to use, default is cyl_cor_single
#' @param verb Depreciated
#' @param ncores If ncores>1, parallelization is used
#' @return Returns a list of cholesky factors of the Vecchia
#' precision matrix corresponding to the years in yearlist
#' @export
build_veccmat_list_grouped<-function(augdata_ordered,grouping_list,
                             yearlist=2007:2016,corrfun=cyl_cor_single,
                             ncores=1,verb=F){
  build_veccmat<-function(year){
    yeardata_ordered=augdata_ordered[augdata_ordered$years==year,]
    groupobj=grouping_list[[year-2007+1]]
    veccmat=create_veccmat_grouped(yeardata_ordered,groupobj,
                                   corrfun=corrfun,ncores=ncores)
    return(veccmat)
  }
  veccmat_list=lapply(yearlist,build_veccmat)
  return(veccmat_list)
}

#'Creat grouped veccmat objects
#'@description Computes cholesky factors of the Vecchia precision matrix corresponding
#'to a specified grouping of conditioning sets. Code inspired by similar code in
#'GPVecchia package, see source code for citation.
#'
#'@param augdata_ordered Ordered data over which to compute the factors. Data
#' should be 'augmented' meaning that each row has the values of the kernel
#' parametmers for that locations; this can be added to a dataframe by calling
#' "augment_data" with a parameter object
#'@param groupobj The grouping information for the Vecchia process
#'@param ncores Number of cores to use, if ncores>1 then parallelization
#'@param corrfun Which function to use for correlations; default cylindrical
#'is used
#'@param covparms Optional, only needed if stationary corrfun is used. For other
#'correlation functions the non-stationary parameters for the kernel convolutions
#'should be included the augdata_ordered dataframe
#'@export
create_veccmat_grouped<-function(augdata_ordered,groupobj,ncores=1,
                                 corrfun=cyl_cor_single,covparms=NULL){
  #this function is based off of the function createU from the GPvecchia
  #package by Katzfuss, Jurek, Zilber, and Gong as found here:
  #https://github.com/katzfuss-group/GPvecchia/blob/master/R/createU.R
  #see paper for full citation
  n=dim(augdata_ordered)[1]

  groupings=groupobj$grouping
  indout=groupobj$indout

  block_entry<-function(NNobj){
    condinds=NNobj[,1]
    inds=sort(unique(c(NNobj)))
    iscond=inds %in% condinds
    if(is.null(covparms)){
      Sigma=outer(inds,inds,corrfun,augdata_ordered)
    }else{
      Sigma=corrfun(covparms,augdata_ordered[inds,])
    }
    A=matrix(0,length(inds),length(condinds))
    A[iscond,]=diag(1,nrow=sum(iscond),ncol=sum(iscond))
    a=solve(chol(Sigma),A)
    a[col(a)<(row(a)-(dim(a)[1]-dim(a)[2]))]=NA #makes easier to remove 0s before sparsity step
    return(a)
  }

  blockentry_list=mclapply(groupings,block_entry,mc.cores=ncores)
  blockentry_vec=unlist(lapply(blockentry_list,function(x) c(x)))
  blockentry_vec=blockentry_vec[!is.na(blockentry_vec)]
  veccmat=Matrix::sparseMatrix(i=indout$rind,
                         j=indout$cind,
                         x=blockentry_vec,
                         dims=c(n,n))
  return(veccmat)
}
