# 08_ohc_trend_posterior_intervals.R
# Calculates posterior credible intervals for yearly ohc values and the
# trend field
# IN: ../MCMC_Output/Output_20k_run/ohc_df.RData
#     ../MCMC_Output/Output_20k_run/mu_slope_df.RData
# OUT: ../FigureScripts/FigureData/posterior_intervals.RData

testmode=T
library(tidyverse)
library(BayesianOHC)

load(file='../MCMC_Output/Output_20k_run/ohc_df.RData',verb=T)
load(file='../MCMC_Output/Output_20k_run/mu_slope_df.RData',verb=T)

#resample ohc values:
burnin_iter=5401
nresamp=100

zalpha=qnorm(.975)
full_resample=F
#faster and should be more accurate to do second rather than
#fully resample
if(full_resample){
  ohc_resampled=ohc_df %>%
    dplyr::filter(iter>burnin_iter)%>%
    mutate(count=nresamp)%>%
    uncount(count)%>%
    mutate(ohc_resamp=rnorm(n(),ohc,sd))
  ohc_intervals_df=ohc_resampled%>%
    group_by(year)%>%
    summarise(ohc_med=median(ohc),
              upper_conf=ohc_med+zalpha*median(sd),
              lower_conf=ohc_med-zalpha*median(sd),
              upper_cred=quantile(ohc_resamp,.975),
              lower_cred=quantile(ohc_resamp,.025))
}else{
  ohc_intervals_df=ohc_df%>%
    group_by(year)%>%
    summarise(ohc_med=median(ohc),
              upper_conf=ohc_med+zalpha*median(sd),
              lower_conf=ohc_med-zalpha*median(sd),
              upper_cred=quantile(ohc+zalpha*sd,.975),
              lower_cred=quantile(ohc-zalpha*sd,.025))
}

yearseq_adj=seq(2006.5,2016.5,length=1000)-2007

#do the same for the trend fields:
trend_intervals_df=mu_slope_df%>%
  expand_grid(yearseq_adj)%>%
  mutate(sd=sqrt(mu0_sd^2+2*corr_term*mu0_sd*slope_sd*yearseq_adj+
                   slope_sd^2*yearseq_adj^2),
         val=mu0+slope*yearseq_adj,
         year=yearseq_adj+2007)%>%
  group_by(year)%>%
  summarise(ymed=median(val),
            upper_slope_conf=ymed+zalpha*median(sd),
            lower_slope_conf=ymed-zalpha*median(sd),
            upper_slope_cred=quantile(val+zalpha*sd,.975),
            lower_slope_cred=quantile(val-zalpha*sd,.025))

if(!test_mode){
  save(ohc_intervals_df,trend_intervals_df,
     file='../MCMC_Output/Output_20k_run/posterior_intervals.RData')
}

