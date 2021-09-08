# 06_extract_samples.R
# Extracts samples from files in the MCMC_Output/Output../Samples
# directory and converts them to matrix form
# IN: Files in the directory '../MCMC_Output/Output_20k_run/Samples/'
# OUT: '../MCMC_Output/Output_20k_run/sample_matrices.RData'
# TEST: Does not overwrite data

test_mode=T
library(dplyr)
library(BayesianOHC)

sampledir=paste('../MCMC_Output/Output_20k_run/Samples/',sep='')
myfiles=list.files(sampledir)
burnin_index=5401

#function for reading in the contents of each file in the Samples directory
get_params_from_file<-function(filename){
  print(filename)
  load(paste(sampledir,filename,sep=''))
  return(param_stores)
}

all_outs=do.call(c,map(myfiles,get_params_from_file))
tts=do.call(c,map(all_outs,'tt'))
end_of_cycle=map(all_outs,'end_of_cycle')
end_of_cycle[sapply(end_of_cycle,is.null)]=F
end_of_cycle=do.call(c,end_of_cycle)

all_params=map(all_outs[end_of_cycle & tts>burnin_index],'stored_params')

#function for converting stored samples into matrix
extract_samples<-function(varname) all_params%>%
  map(paste('basis',ifelse(varname%in%c('slope','mu0'),'mu','fields'),sep='_'))%>%
  map(paste('basis',varname,sep='_'))%>%
  do.call(cbind,.)

sample_matrices=map(default_varname_list,extract_samples)%>%
  set_names(default_varname_list)

#store mean hyper-mean as it is the only hyperparameter that is updated:
sample_matrices$hyper_mu0_mean=all_params%>%
  map(~.$hyperparam_list$mu0$mu)%>%
  do.call(c,.)

if(!test_mode){
  save(sample_matrices,
       file='../MCMC_Output/Output_20k_run/sample_matrices.RData')
}
