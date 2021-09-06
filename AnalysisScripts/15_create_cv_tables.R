# 15_create_cv_tables.r
# Creates LOFO and LOWO cross validation tables corresponding to Table 3 and
# Table S2 respectively. Variable cv_type can be "lofo" or "lowo"
# IN: For each model_name, files of the name
#     ../ValidationData/lofodf_<model_name>.RData
#     (and analogous for LOWO)
# OUT: prints xtable object
# TEST: N/A, does not save output

library(tidyverse)
library(BayesianOHC)
library(xtable)
library(verification)

model_names=c('fullmodel','stat_phi','stat_theta_lat','stat_theta_lon',
              'stat_nugget','stat_all',
              'levitus','chordal','spacetime','chordal_nonstat')%>%
  set_names(.)

cv_type='lofo'
filenames=map(model_names,~paste('../ValidationData/',cv_type,'df_',.,'.RData',sep=''))
cv_dfs=map(filenames,~get(load(.)))

texnames=list('fullmodel'='Fully Nonstationary',
              'stat_nugget'='Stationary $\\sigma^2/\\phi$',
              'stat_theta_lat'='Stationary $\\theta_{\\textnormal{lat}}$',
              'stat_theta_lon'='Stationary $\\theta_{\\textnormal{lon}}$',
              'stat_phi'='Stationary $\\phi$',
              'chordal'='Stationary Chordal',
              'spacetime'='Stationary Spatiotemporal',
              'levitus'='Levitus et al.',
              'stat_all'='Fully Stationary',
              'chordal_nonstat'='Nonstationary Chordal')

program=sapply(model_names,function(x) 'BayesOHC')
program['chordal']='GpGp'
program['spacetime']='GpGp'

myrmses=sapply(cv_dfs,get_rmse)
ref_rmse=get_rmse(cv_dfs$levitus)
formatfun3=function(x) formatC(x,digits=3,format='f',flag='#')
formatfun2=function(x) formatC(x,digits=2,format='f',flag='#')
cv_table=data.frame(program=program[model_names],
              rmse=myrmses,
             crps=sapply(cv_dfs,get_crps))%>%
  mutate(improvement=paste(formatfun2(100*(1-rmse/ref_rmse)),'\\%',sep=''))%>%
  arrange(rmse)%>%
  mutate(rmse_display=formatfun3(rmse),
         crps=formatfun3(crps))%>%
  dplyr::select(!c(rmse))%>%
  relocate(rmse_display,.after=program)%>%
  relocate(crps,.after=improvement)%>%
  set_names(c('Program','RMSE (GJ/m$^2$)','Improvement/Ref','CRPS'))%>%
  `rownames<-`(texnames[rownames(.)])
cv_table
print(xtable(cv_table),sanitize.text.function=function(x){x})

