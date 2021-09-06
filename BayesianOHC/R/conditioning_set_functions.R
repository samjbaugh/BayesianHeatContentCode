#' Computes list of grouped conditioning sets for each year in the given data
#'
#' @param locations Ordered data on which to compute nearest neighbors
#' @param m Number of nearest neighbors to use in constructing the ungrouped
#' vecchia conditioning sets.
#' @param yearlist List of years over which to compute the group sets.
#' @return Returns a list of length yearlist where each entry contains the
#' output of compute_grouped_conditioning_sets on the data for the
#' corresponding year.
#' @param verb T/F, if T displays a progress bar
#' @param ncores Number of cores to use, if ncores>1 parallelization is
#' @importFrom parallel mclapply
#' @export
compute_grouping_list<-function(locations,m,yearlist=2007:2016,
                                verb=F,ncores=1){
  mylist=list()
  for(year in yearlist){
    yeardata_ordered=locations[locations$years==year,]
    myNNarray=compute_nnarray(yeardata_ordered,m=m,verb=verb,ncores=ncores)
    mygroupobj=compute_grouped_conditioning_sets(myNNarray)
    mylist[[year-min(yearlist)+1]]=mygroupobj
  }
  return(mylist)
}

#' Computes array of nearest neighbors to be used for Vecchia process
#' conditioning sets
#'
#' @param locations Ordered data on which to compute nearest neighbors
#' @param m Number of nearest neighbors
#' @param verb T/F, if T displays a progress bar
#' @param ncores Number of cores to use, if ncores>1 parallelization is
#' used
#' @return Returns a matrix of size (nlocs,m) where nlocs is the number of
#' locations specified in "ordered_input".
#' @export
compute_nnarray<-function(locations,m,verb=F,ncores=1)
{
  if(is.null(locations$obs)){
    locations$obs=T
  }
  obs_df=locations[locations$obs,]
  myn=dim(locations)[1]
  obs_df=locations

  nobs=dim(obs_df)[1]
  NNarray <- matrix(NA,myn,m+1)
  nnarray_entry<-function(j){
    if(j==1){return(c(1,rep(NA,m)))}
    jind=min(j-1,nobs)
    distsamps=1:jind
    distvec=c(gen_pdist_mat_cylindrical(obs_df[distsamps,c('lat_rad','lon_rad'),drop=FALSE],
                                           locations[j,c('lat_rad','lon_rad'),drop=FALSE]) )

    return(c(j,sort(c(order(distvec)[1:min(m,(j-1))],
                      rep(NA,max(0,(m)-(j-1)))),decreasing=TRUE,na.last=T)))
  }
  NNarray=do.call(rbind,mclapply(1:myn,nnarray_entry,mc.cores=ncores))
  return(NNarray)
}


#' Computes grouped conditioning sets for use in the Vecchia approximation.
#'
#' @param NNarray Matrix containing indices of nearest neighbors; generally
#' returned from compute_nnarray
#' @param trivial If true returns grouping object with no grouping done
#' @export
compute_grouped_conditioning_sets<-function(NNarray,trivial=F){
  if(trivial){
    grouping_unord2=lapply(1:dim(NNarray)[1],function(x) NNarray[x,,drop=FALSE])
  }else{
    grouping_unord=compute_nn_clusters(NNarray)
  }
  grouping=lapply(grouping_unord,function(x) x[order(x[,1]),,drop=FALSE])

  frc=function(NNobj){
    condinds=NNobj[,1]
    inds=sort(unique(c(NNobj)))
    ninds=length(inds)
    nconds=length(condinds)
    cind=unlist(sapply(1:nconds,function(jj) rep(condinds[jj],ninds-(nconds-jj))))
    rind=unlist(sapply(1:nconds,function(jj) inds[1:(ninds-(nconds-jj))]))
    return(data.frame('cind'=cind,'rind'=rind))
  }

  indout=do.call(rbind,lapply(grouping,frc))

  return(list(grouping=grouping,indout=indout))
}

#' Computes grouping object for full Cholesky
#'
#' @param nobs Number of observations
#' @export
groupinfo_full_chol<-function(nobs){
  grouping=list(t(sapply(1:nobs,function(j) c(j:1,rep(NA,nobs-j+1)))))
  frc=function(NNobj){
    condinds=NNobj[,1]
    inds=sort(unique(c(NNobj)))
    ninds=length(inds)
    nconds=length(condinds)
    cind=unlist(sapply(1:nconds,function(jj) rep(condinds[jj],ninds-(nconds-jj))))
    rind=unlist(sapply(1:nconds,function(jj) inds[1:(ninds-(nconds-jj))]))
    return(data.frame('cind'=cind,'rind'=rind))
  }
  indout=do.call(rbind,lapply(grouping,frc))
  return(list(grouping=grouping,indout=indout))
}


#' Finds clusters of sets of nearest neighbors for grouping sets. Function
#' used internally and is adapted from GpGP package (see in-function comments
#' for citation details).
#'
#' @param NNarray Matrix continaing indices of nearest neighbors; generally
#' returned from compute_nnarray
#' @return Returns a list of clustered arrays of nearest neighbors of length
#' m, where m is the second dimension of NNarray
compute_nn_clusters<-function(NNarray){
    ## Citation: The following code segment in brackets is adapted from the
    ## source code of the function "group_obs" of the GpGp package version 0.4.0,
    ## accessed from https://github.com/cran/GpGp/blob/master/R/nearest_neighbor_functions.R
    ## authorized under the MIT + license. Credit goes to the authors
    ## Joseph Guinness, Matthias Katzfuss, and Fahmy Youssef
    n = nrow(NNarray)
    m = ncol(NNarray)-1

    clust = vector("list",n)
    for(j in 1:n) clust[[j]] <- j
    for( ell in 2:(m+1) ){  # 2:(m+1)?
      sv = which( NNarray[,1] - NNarray[,ell] < n )
      for(j in sv){
        k = NNarray[j,ell]
        if( length(clust[[k]]) > 0){
          nunique = length(unique(c(NNarray[c(clust[[j]],clust[[k]]),])))
          if( nunique^2 <= length(unique(c(NNarray[clust[[j]],])))^2 +
              length(unique(c(NNarray[clust[[k]],])))^2 ) {
            clust[[j]] = c(clust[[j]],clust[[k]])
            clust[[k]] = numeric(0)
          }
        }
      }
    }
    zeroinds <- unlist(lapply(clust,length)==0)
    clust[zeroinds] <- NULL
    retval= lapply(clust,function(inds) NNarray[inds,,drop=FALSE])
    return(retval)
}
