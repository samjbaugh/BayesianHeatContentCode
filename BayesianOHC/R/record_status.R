#' Records status of mcmc sampler
#' @description Records the status of the mcmc sampler through appending to the file
#' specified in "status_filename".
#'
#' @param status_filename Location to write the sampler status.
#' @param iteration The current iteration.
#' @param current_likelihood The current likelihood.
#' @param current_prior The current prior value.
#' @param accepts The current vector of Metropolis-Hastings acceptances.
#' @param variance_scaling_factor The current vector of scaled proposal
#' variances.
#' @param start_time The original start time for the sampler
#' @return Returns a data frame that is a copy of the input data but with new
#' parameter values for each of the fields specified in varname_list
#' @export
record_status<-function(status_filename,iteration,current_likelihood,
                        current_prior,accepts,variance_scaling_factor,
                        start_time)
{
  sink(file=status_filename,append=T)
  cat(sprintf('\n--------------------------------------------------------\n'))
  print(data.frame('iteration'=iteration,
                   'current_likelihood'=current_likelihood,
                   'current_prior'=current_prior,
                   'time'=Sys.time(),
                   'timediff'=Sys.time()-start_time))
  toprint=apply(accepts[1:iteration,],2,mean)
  names(toprint)=sapply(names(accepts),function(x) paste('accept_',x,sep=''))
  print(toprint)
  print(data.frame(variance_scaling_factor))
  cat(sprintf('--------------------------------------------------------'))
  sink()
}
