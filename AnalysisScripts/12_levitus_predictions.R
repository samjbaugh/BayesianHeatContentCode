# 12_levitus_predictions.R
# Uses the method described in Levitus et al. (2012) to produce ocean heat 
# content anomalies to be used in cross validation, as well as predictions
# of ocean heat content by year. Procedure is described in detail in the
# comments below. This scripts also compares trend posteriors with the
# full BayesianOHC model fit.
# IN: ../MCMC_Input/argo_data_january.RData
#     ../MCMC_Output/Output_20k_run/trend_resamp_intonly.RData
# OUT: ../ValidationData/levitus_anomalies.RData
#      ../ValidationData/levitus_ohc.RData
#      Also prints an xtable object corresponding to Table S1
# TEST: Does not overwrite output

# Procedure:
# Originally described in the Levitus et al. (2012) supplementary material found here:
# https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2012GL051106&file=grl29030-sup-0003-txts01.pdf
# In summary we do the following:
#1. Identify points within an influence region around the desired point
#2. Each location is given a "first guess" of zero anomaly, which is corrected
#   by adding   weighted sum of observations within the influence region (see
#   equation 3 in the above linked document). The weights are w=exp(-Er^2/R^2)
#   where r is the distance between the points in km, R is the radius of the
#   influence region (set to be 666km), and E=4.
#3. Standard errors are calculated with the following formula (see equation
#   (7) in the above linked document):
#     sigmaA=sigma0*sqrt(sum((w_i/W)^2)+2*sum(sum(wiwj))) where sigma0 is the
#     empirical standard deviation of the "corrected" values
#     within the influence region
testmode=T
library(tidyverse)
library(BayesianOHC)
library(geodist)

load('../MCMC_Input/argo_data_january.RData',verb=T)

#Compute anomalies; this is equivalent to the "first pass" described in 
# the Levitus et al. supplementayr materials
roundfun=function(x) round(x-1/2)+1/2
firstpass_obsdf=bind_rows(argo_data_ordered%>%mutate(obs=T),
                      gen_masked_grid(1,1)%>%mutate(obs=F))%>%
  mutate(mu0=0,sd=sd(argo_data_ordered$vhc_obs))
gridmeans=levitus_preds(firstpass_obsdf$obs==F,obspred_df=firstpass_obsdf,verb=T)%>%
  mutate(lon_round=roundfun(lon_degrees),
         lat_round=roundfun(lat_degrees), #round((lat_degrees-1/2)/2)*2+1/2,
                               mu0=validval,sd=validsd)

#add computed anomalies to dataframe
argo_data_levitus_anoms=argo_data_ordered%>%
  mutate(lon_round=roundfun(lon_degrees),
         lat_round=roundfun(lat_degrees))%>%
  left_join(gridmeans%>%dplyr::select('lat_round','lon_round','mu0','sd'),
            c('lat_round','lon_round'))%>%
  mutate(mu0=ifelse(is.na(mu0),mean(vhc_obs),mu0))%>%
  mutate(sd=ifelse(is.na(sd),mean(vhc_obs),sd))%>%
  mutate(z=vhc_obs-mu0)

if(!testmode){
  save(argo_data_levitus_anoms,gridmeans,
       file='../ValidationData/levitus_anomalies.RData')
}

###Compute OHC predictions and uncertainties
obspred_df=bind_rows(argo_data_ordered%>%mutate(obs=T),
                     gen_masked_grid(2,2)%>%mutate(obs=F,years=NA))%>%
  mutate(lon_round=roundfun(lon_degrees),
         lat_round=roundfun(lat_degrees))%>%
  left_join(gridmeans%>%dplyr::select('lat_round','lon_round','mu0','sd'),
            c('lat_round','lon_round'))%>%
  dplyr::filter(!is.na(mu0))%>%
  mutate(mu0=ifelse(is.na(mu0),mean(vhc_obs,na.rm=T),mu0))%>%
  mutate(sd=ifelse(is.na(sd),mean(vhc_obs),sd))%>%
  mutate(z=vhc_obs-mu0)

ohcdf_list=list()
for(myyear in 2007:2016){
  print(myyear)
  year_obspred_df=subset(obspred_df,years==myyear | is.na(years))
  preds_out=levitus_preds(pred_indices=!year_obspred_df$obs,
                          obspred_df=year_obspred_df,distfun=geodist)
  ohcdf_list[[myyear-2007+1]]=preds_out
}

scalarvec=gen_deltas_from_grid(ohcdf_list[[myyear-2006]],res=deg2rad(2))
levitus_ohc=data.frame(year=yearlist,
                       ohc=sapply(ohcdf_list,function(x) sum(x$validval*scalarvec)),
                       ohc_sd=sapply(ohcdf_list,function(x) sqrt(sum(x$validsd^2*scalarvec^2))))
if(!test_mode){
  save(levitus_ohc,file='../ValidationData/levitus_ohc.RData')
}

#Create Table S1; load levitus_ohc and trend_resamp_intonly to compare
#posteriors between the two methods:
load('../ValidationData/levitus_ohc.RData',verb=T)
load('../MCMC_Output/Output_20k_run/trend_resamp_intonly.RData',verb=T)
levitus_lm=lm(data=levitus_ohc,ohc~year)
levitus_lm$coefficients[2]
confint(levitus_lm)[2,]
levitus_lm_summary=summary(levitus_lm)
quantile(resamp_int_trends,c(0.025,.5,.975))
#by the way calculate p-value of BayesianOHC posterior:
mean(resamp_int_trends<0)

formatfun=function(x) {
  x2=as.numeric(x/1e12)
  paste(formatC(x2, format = "f", digits = 3),'$\\times 10^{21}$',sep='')
}

levitusTable=data.frame(low_bound=formatfun(c(confint(levitus_lm)[2,1],quantile(resamp_int_trends,0.025))),
           est=formatfun(c(levitus_lm$coefficients[2],median(resamp_int_trends))),
           upper_bound=formatfun(c(confint(levitus_lm)[2,2],quantile(resamp_int_trends,.975))),
           width=formatfun(c(diff(confint(levitus_lm)[2,]),diff(quantile(resamp_int_trends,c(0.025,.975))))),
           pvalue=formatC(c(levitus_lm_summary$coefficients[2,4],mean(resamp_int_trends<0)),format='f',digits=5))
names(levitusTable)=c('2.5\\%','Median','97.5\\%','Width','p-value')
rownames(levitusTable)=c('Levitus et al.','BayesianOHC')

library(xtable)
print(xtable(levitusTable),sanitize.text.function=function(x){x})
