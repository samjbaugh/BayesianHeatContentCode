#' Runs MCMC sampler for the cylindrical model with linear trend term
#'
#' @param ordered_data Data for running the sampler; should be already ordered.
#' @param initparams Parameters for initializing the sampler.
#' @param M number of iterations to run the sampler.
#' @param m Can either be the string "chol" indicating that the full
#' Cholesky should be used, or an integer indicating the number of Vecchia
#' neighbors that should be used.
#' @param yearlist The list of years to be used.
#' @param pred_locs Locations to compute predictions
#' @param grouping_list List of observation location groupings for Vecchia
#' approximation, if not included the list will be computed using the value of m.
#' @param grouping_list_preddf List of groupings of prediction locations
#' (as specified in pred_locs); if empty this will be computed in the function
#' @param outdir Path to directory for storing output.
#' @param pred_iter If do_predictions=T the sampler will compute predictions
#' for every pred_iter iterations
#' @param varname_list List of variable names to sample.
#' @param corrfun Correlation fun to use, defaults to cylindrical
#' @param var_sample_order Order of variables names for sampling. Note that mu0
#' and slope should be specified as "mu" in this list as they are sampled together
#' from their joint marginal posterior distribution
#' @param stationary_varnames List containing which of the specified variables
#' should be maintained as stationary.
#' @param region_name Name of the region to restrict to; see documentation
#' of "in_region" for available  options. Data does not need to be restricted
#' to the region beforehand.
#' @param continue_from_previous T/F, if T continue from previous iteration.
#' @param run_label Label used in saving the output. If label is reused and
#' continue_from_previous=F then the previous run will be overwritten
#' @param save_iter Every save_iter iterations, the samples since the last
#' checkpoint will be saved in the output directory. Default is to not save
#' samples and rather return all of the samples at the end.
#' @param ncores Number of cores to be used; if ncores>1 then parallelization
#' will be used
#' @param linkfuns List of link functions for the variables specified
#' in varname_list.
#' @param invlinkfuns Inverse link functions for the variables specified
#' in varname list.
#' @param ret_last_veccmat T/F, if T returns the last veccmat list
#' @param plot_iter At each plot_iter iterations figures showing the
#' parameter field maps will be saved in the output directory (note: plotting
#' not currently implemented)
#' @return Returns the last M-save_iter samples, or all M samples if
#' save_iter>M. If save_iter<M, earlier samples are stored every save_iter
#' iterations in the directory yspecified by outdir
#' @import GPvecchia
#' @importFrom parallel mclapply
#' @export
run_mcmc_sampler<-function(ordered_data,initparams,
                           M=20000,m="chol",yearlist=2007:2016,
                           pred_locs=NULL,grouping_list=NULL,
                           grouping_list_preddf=NULL,outdir='../MCMC_Output/',
                           varname_list=NULL,corrfun=cyl_cor_single,
                           var_sample_order=NULL,
                           stationary_varnames=NULL,
                           region_name='globe',continue_from_previous=F,
                           run_label='mcmc',
                           linkfuns=default_linkfuns,ret_last_veccmat=F,
                           invlinkfuns=default_invlinkfuns,
                           pred_iter=Inf,save_iter=Inf,
                           plot_iter=Inf,ncores=1){
  ############Initialize Sampler Settings ###############
  #not currently implemented:
  plot_figures=F #plot_iter<M

  if(is.null(varname_list)){
    varname_list=names(initparams$hyperparam_list)
  }
  basis_varname_list=sapply(varname_list,function(x) paste('basis',x,sep='_'))
  if(is.null(var_sample_order)){
    var_sample_order=varname_list[varname_list!='slope']
  }
  if(is.null(stationary_varnames)){
    stationary_varnames=names(initparams$stationary_params)
  }

  ordered_data=in_region(ordered_data,region_name)

  #Initialize output directories#

  if(is.infinite(save_iter)){
    #If save not specified, assume output should be suppressed
    storage_dirname='/dev/null'
    status_filename='/dev/null'
  }else{
    if(!file.exists(outdir)) {dir.create(file.path(outdir))}
    storage_dirname=paste(outdir,'/Output_',run_label,sep='')
    samples_dirname=paste(outdir,'/Output_',run_label,'/Samples/',sep='')
    if(!file.exists(storage_dirname)) {dir.create(file.path(storage_dirname))}
    if(!file.exists(samples_dirname)) {dir.create(file.path(samples_dirname))}
    status_filename=paste(storage_dirname,'/output_log.txt',sep='')
  }

  if(plot_figures){
    figure_dirname=paste(storage_dirname,'/Figures/',sep='')
    if(!file.exists(figure_dirname)) {dir.create(file.path(figure_dirname))}
    if(!file.exists(paste(figure_dirname,'/parameter_fields',sep=''))) {dir.create(file.path(figure_dirname,'parameter_fields'))}
    if(!file.exists(paste(figure_dirname,'/validation_images',sep=''))) {dir.create(file.path(figure_dirname,'validation_images'))}
    if(!file.exists(paste(figure_dirname,'/krigged_values',sep=''))) {dir.create(file.path(figure_dirname,'krigged_values'))}
    if(!file.exists(paste(figure_dirname,'/mus',sep=''))) {dir.create(file.path(figure_dirname,'mus'))}
  }

  do_predictions=!is.null(pred_locs)
  if(do_predictions){
    pred_locs=in_region(pred_locs,region_name)
    scalarvec=gen_deltas_from_grid(pred_locs)
  }

  prop_batch_size=100
  mu_indices=which(varname_list%in%c('mu0','slope'))

  if(is.null(grouping_list)){
    if(m=="chol"){
    nyears=sapply(yearlist,function(year) sum(ordered_data$years==year))
    grouping_list=lapply(nyears,groupinfo_full_chol)
    }else{
    grouping_list=compute_grouping_list(ordered_data,m=m)
    }
  }
  if(do_predictions & is.null(grouping_list_preddf)){
    pred_df=ordered_data[,c(names(pred_locs),'years')]
    pred_df$obs=T
    pred_locs$obs=F
    for (year in 2007:2016) {
      pred_locs$years=year
      pred_df=rbind(pred_df,pred_locs)
    }
    grouping_list_preddf=compute_grouping_list(pred_df,m=m)
  }

  nonstationary_varnames=varname_list[!(varname_list%in%stationary_varnames)]
  num_steps_per_iter=length(var_sample_order)
  Mii=(M)*num_steps_per_iter+1
  nonstat_indices=which(!(varname_list%in%stationary_varnames |
                            varname_list%in%c('mu0','slope')))
  calc_prior<-function(myparams){
    return(-sum(unlist(c(myparams$basis_fields[,basis_varname_list[nonstat_indices]]^2,0)))/2-
             sum(unlist(c(myparams$basis_mu[,basis_varname_list[mu_indices]]^2,0)))/2-
             ifelse(length(stationary_varnames)>0,
                    sum(unlist(c(sapply(stationary_varnames,
                                        function(vv) (invlinkfuns[[vv]](myparams$stationary_params[[vv]])-
                                                        myparams$hyperparam_list[[vv]]$mu)^2/
                                          (2*sqrt(myparams$hyperparam_list[[vv]]$phi)))),0)),0))
  }

  ############Load Initial Configuration###############
  if(!continue_from_previous)
  {
    current_params=get_region_params(initparams,region_name,
                                     varname_list,invlinkfuns)
    if(length(stationary_varnames)==0){
      current_params$stationary_params=list()
    }
    current_augdata=augment_data_parallel(ordered_data,current_params,
                                          varname_list,ncores=ncores)
    current_veccmat_list=build_veccmat_list_grouped(current_augdata,grouping_list,
                                                    corrfun=corrfun,ncores=ncores)

    if(region_name!='globe'){
      #update mean params to posterior mean; only needs to be done here if
      #restricting mean field to a sub-region of the globe
      current_params=sample_mean_trend(current_augdata,
                                       current_params,return_posterior_mean=T,
                                       current_veccmat_list,return_full_posterior=F,
                                       ncores=ncores)
      current_augdata=augment_data(current_augdata,current_params,
                                    c('mu0','slope'),linkfuns)
    }

    current_likelihood=sum(vecc_lik_over_years(current_augdata,current_veccmat_list))
    current_prior=calc_prior(current_params)

    sample_list=list(list(stored_params=current_params,
                          current_likelihood=current_likelihood,
                          current_prior=current_prior,
                          ii=1,
                          tt=0,
                          myvarname=NA,
                          end_of_cycle=0))

    pred_stores=list()
    value_stores=data.frame(current_likelihoods=rep(NA,Mii),
                            current_priors=rep(NA,Mii),
                            proposed_likelihoods=rep(NA,Mii),
                            proposed_priors=rep(NA,Mii),
                            is_end_of_cycle=rep(NA,Mii),
                            ii=rep(NA,Mii),
                            tt=rep(NA,Mii))
    mu_slope_df=data.frame(mu0_mean=rep(NA,M),
                           mu0_val=rep(NA,M),
                           slope_mean=rep(NA,M),
                           slope_var=rep(NA,M),
                           covar=rep(NA,M),
                           tt=rep(NA,M))
    value_stores[1,]=c(current_likelihood,current_prior,
                       current_likelihood,current_prior,TRUE,1,1)

    accepts=do.call(cbind,lapply(basis_varname_list[-mu_indices],function(vv){
      a=data.frame(c(0,rep(NA,max(M,1)-1)))
      names(a)=vv
      return(a)
    }))
    starti=2

    variance_scaling_factor=as.list(rep(2,4))
    names(variance_scaling_factor)=basis_varname_list[-mu_indices]
    nbasis=dim(current_params$basis_fields)[1]

    start_time=Sys.time()
    record_status(status_filename,1,
                  current_likelihood,current_prior,
                  accepts,variance_scaling_factor,start_time)
  }else
  {
    if(!file.exists(paste(storage_dirname,'/current_values.RData',sep=''))){
      print('Error: cannot continue from previous (invalid directory)')
      return(NA)
    }
    load(paste(storage_dirname,'current_values.RData',sep=''),verbose=T)
    value_stores$ii=NA
    value_stores$ii[1:ii]=1:ii
    value_stores$tt=value_stores$ii%/%num_steps_per_iter
    if(!exists("mu_slope_df")){
      mu_slope_df=data.frame(mu0_mean=rep(NA,M),
                             mu0_val=rep(NA,M),
                             slope_mean=rep(NA,M),
                             slope_var=rep(NA,M),
                             covar=rep(NA,M))
    }

    sample_list=list(list(stored_params=current_params,
                           ii=ii,
                           tt=tt))
    starti=ii+1 #start at next iternum
    print(paste('Existing output found, resuming at iteration',tt))

    current_augdata=augment_data_parallel(ordered_data,current_params,
                                         varname_list,linkfuns,ncores=ncores)
    for(varname in stationary_varnames){
      current_augdata[[varname]]=current_params$stationary_params[[varname]]
    }
    current_likelihood=sum(vecc_lik_over_years(current_augdata,current_veccmat_list))
    current_prior=calc_prior(current_params)

    nbasis=dim(current_params$basis_fields)[1]
    start_time=Sys.time()

  }
  ############Begin Sampler Loop###############

  ##steps per iter guide (8 steps):
  # ii%%10=1, update basis_phi
  # ii%%10=2, update basis_nugget
  # ii%%10=3, update basis_theta_lat
  # ii%%10=4, update basis_theta_lon
  # ii%%10=5, update basis_mu0 and basis_slope simultaneously
  for(ii in starti:Mii)
  {
    if(Mii<starti){break}

    tt=((ii-2)%/% num_steps_per_iter) + 1
    step_id=(ii-2)%%num_steps_per_iter + 1

    end_of_cycle=(step_id==(num_steps_per_iter))
    varname=var_sample_order[step_id]

    if(end_of_cycle)
    {
      print(tt)
    }

    if(varname %in% c('theta','theta_lat','theta_lon','nugget','phi')) #update parameter fields
    {
      basis_varname=paste('basis',varname,sep='_')
      if(varname %in% stationary_varnames){
        proposed_params=current_params
        proposed_params$stationary_params[[varname]]=
          linkfuns[[varname]](invlinkfuns[[varname]](current_params$stationary_params[[varname]])+
          rnorm(1,sd=variance_scaling_factor[[basis_varname]]))
        proposed_augdata=current_augdata
        proposed_augdata[[varname]]=proposed_params$stationary_params[[varname]]

        if(varname!='phi'){ #don't need to recreate veccmat for phi
          proposed_veccmat_list=tryCatch(
            build_veccmat_list_grouped(proposed_augdata,grouping_list,
                                       corrfun=corrfun,ncores=ncores),
                                         error=function(e) return(NULL))
        }else{
          proposed_veccmat_list=current_veccmat_list
        }
        if(is.null(proposed_veccmat_list)){
          proposed_likelihood=-Inf
          }else{
          proposed_likelihood=sum(vecc_lik_over_years(proposed_augdata,proposed_veccmat_list))
        }
        proposed_prior=calc_prior(proposed_params)

        alpha=exp(proposed_likelihood-current_likelihood+proposed_prior-current_prior)
        randnum=runif(1) #one for each k
        accepted=randnum<alpha
        accepts[[basis_varname]][tt]=accepted
        if(accepted)
        {
          current_params=proposed_params
          current_likelihood=proposed_likelihood
          current_augdata=proposed_augdata
          current_prior=proposed_prior
          current_veccmat_list=proposed_veccmat_list
        }

      }else{
        #construct proposal
        proposed_params=current_params
        proposed_params$basis_fields[[basis_varname]]=
          current_params$basis_fields[[basis_varname]]+
          rnorm(nbasis,sd=variance_scaling_factor[[basis_varname]])
        proposed_augdata=augment_data(current_augdata,proposed_params,
                                      varname,linkfuns)

        if(varname!='phi'){ #don't need to recreate veccmat for phi
          proposed_veccmat_list=
            tryCatch(build_veccmat_list_grouped(proposed_augdata,grouping_list,
                                                corrfun=corrfun,ncores=ncores),
                                         error=function(e) return(NULL))
        }else{
          proposed_veccmat_list=current_veccmat_list
        }
        if(is.null(proposed_veccmat_list)){
          proposed_likelihood=-Inf
        }else{
          proposed_likelihood=sum(vecc_lik_over_years(proposed_augdata,proposed_veccmat_list))
        }
        proposed_prior=calc_prior(proposed_params)

        #run Metropolis-Hastings step
        alpha=exp(proposed_likelihood-current_likelihood+proposed_prior-current_prior)
        randnum=runif(1) #one for each k
        accepted=randnum<alpha
        accepts[[basis_varname]][tt]=accepted
        if(accepted)
        {
          current_params=proposed_params
          current_likelihood=proposed_likelihood
          current_augdata=proposed_augdata
          current_prior=proposed_prior
          current_veccmat_list=proposed_veccmat_list
        }
        #update variance scaling factor based on last batch_len iterations
      if(tt>4)
      {
        ideal_acceptance_rate=.23
        delta=
          ifelse(mean(accepts[[basis_varname]][max(1,tt-prop_batch_size):tt])<
                   ideal_acceptance_rate,-min(0.1,1/sqrt(tt)),min(0.1,1/sqrt(tt)))
        variance_scaling_factor[[basis_varname]]=
          variance_scaling_factor[[basis_varname]]*exp(delta)
        if(variance_scaling_factor[[basis_varname]]>2){
          variance_scaling_factor[[basis_varname]]=2
        }
      }
    }
    }else if(varname%in%c('mu0','slope'))
    {
      sample_mean_out=sample_mean_trend(current_augdata,current_params,
                                        current_veccmat_list,
                                        return_full_posterior=T,ncores=ncores)
      proposed_params=sample_mean_out$retparams

      proposed_augdata=augment_data(current_augdata,proposed_params,
                                    c('mu0','slope'),linkfuns)

      proposed_likelihood=sum(vecc_lik_over_years(proposed_augdata,current_veccmat_list))
      proposed_prior=calc_prior(proposed_params)

      #always accept since its the exact Gibbs distribution
      current_params=proposed_params
      current_likelihood=proposed_likelihood
      current_augdata=proposed_augdata
      current_prior=proposed_prior

      if(do_predictions){
        muslope_dists1=get_mu_integrated_posterior(sample_mean_out,proposed_params,
                                                   pred_locs=pred_locs)
        mu_slope_df[tt,]=data.frame(muslope_dists1[1:5],iter=tt)
      }
    }

    #####Storage, Logs, and Plots####
    #store values and parameters first. Record status if end of cycle
    value_stores[ii,]=data.frame(current_likelihoods=current_likelihood,
                                 current_priors=current_prior,
                                 proposed_likelihoods=proposed_likelihood,
                                 proposed_priors=proposed_prior,
                                 is_end_of_cycle=end_of_cycle,
                                 ii=ii,
                                 tt=tt)
    storage_i=length(sample_list)
    sample_list[[storage_i+1]]=list(stored_params=current_params,
                                     current_likelihood=current_likelihood,
                                     current_prior=current_prior,
                                     ii=ii,
                                     tt=tt,
                                     myvarname=varname,
                                     end_of_cycle=end_of_cycle)
    if(end_of_cycle){
      if(!is.infinite(save_iter)){
        print(paste('Recording'))
        record_status(status_filename,tt,
                      current_likelihood,current_prior,
                      accepts,variance_scaling_factor,start_time)
      }
    }

    if(do_predictions & end_of_cycle & (tt%%pred_iter==0)){
      print('Computing Predictions:')
      pred_ret=compute_lincomb_predictions(obsdata=current_augdata,
                                          pred_locs=pred_locs,
                                          myparams=current_params,
                                          grouping_list_preddf=grouping_list_preddf,
                                        scalarvec = scalarvec,
                                        ncores=ncores)
      pred_ret$iter=tt
      pred_stores[[tt]]=pred_ret
      #rename if not predicting ohc?
      ohc_df=do.call(rbind,pred_stores)
      save(ohc_df,file=paste(storage_dirname,'/ohc_df.RData',sep=''))
    }
    if(end_of_cycle & (tt%%save_iter==0))
    {
      print(paste('Saving'))
      #This is for if you want to restart the sampler where it left off:
      save(current_params,
           current_veccmat_list,
           current_augdata,
           value_stores,
           accepts,
           variance_scaling_factor,
           tt,ii,save_iter,num_steps_per_iter,
           # mean_field,mean_hyper,C_hyper,C0_hyper,
           file=paste(storage_dirname,'/current_values.RData',sep=''))
      save(mu_slope_df,file=paste(storage_dirname,'/mu_slope_df.RData',sep=''))
      save(sample_list,
           file=paste(samples_dirname,'/stored_values_iter_',tt,'.RData',sep=''))
      rm(sample_list) #free up memory
      sample_list=list()
    }
    if(end_of_cycle & (tt%%plot_iter==0))
    {
      #plotting not currently implemented
      print(paste('plotting not currently implemented'))
      # plot_summary_maps_at_obslocs(current_augdata,
      #                              figure_dirname=figure_dirname,
      #                              iternum=tt,region_name=region_name)
      # nlik=sum(!is.na(value_stores$current_likelihoods))
      # likplot=ggplot()+geom_line(aes(x=(1:nlik),
      # y=value_stores$current_likelihoods[1:nlik]))
      # ggsave(paste(out_dir,'/likplot.png',sep=''),likplot)


    }
  }
  if(ret_last_veccmat){
    return(list('sample_list'=sample_list,
                'last_veccmat'=current_veccmat_list))
  }else{
    return(sample_list)
  }
}
