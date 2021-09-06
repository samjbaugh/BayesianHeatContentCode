testmode=T
library(BayesianOHC)
figure_directory='../Figures'
library(tidyverse)

load('../MCMC_Output/Output_20k_run/posterior_intervals.RData',verb=T)

pcredible=ggplot(data=ohc_intervals_df)+
  geom_line(aes(x=year,y=ohc_med))+
  geom_line(data=trend_intervals_df,aes(x=year,y=ymed),linetype='dotted',col='black')+
  geom_ribbon(data=trend_intervals_df,
              aes(x=year,ymin=lower_slope_conf,ymax=upper_slope_conf,col='Conf Int'),
              alpha=.1)+
  geom_line(data=trend_intervals_df,aes(x=year,y=lower_slope_cred,col='Cred Int'),
            alpha=1,linetype='dashed')+
  geom_line(data=trend_intervals_df,aes(x=year,y=upper_slope_cred,col='Cred Int'),
            alpha=1,linetype='dashed')+
  geom_errorbar(aes(x=year-.15,ymin=lower_conf,ymax=upper_conf,col='Conf Int'),
                stat="identity",alpha=.9,width=.4)+
  geom_errorbar(aes(x=year+.15,ymin=lower_cred,ymax=upper_cred,col='Cred Int'),
                stat="identity",alpha=.9,width=.4)+
  scale_y_continuous(name=bquote('OHC (10'^25*' Joules)'),
                     breaks=seq(1513e22,1531e22,by=4e22)/1e9,
                     labels=seq(1.513,1.531,by=.004))+
  scale_x_continuous(breaks=seq(2007,2016,by=1),expand=c(0,0),name='Year')+
  scale_color_manual(values=c('red','blue'))+
  theme(legend.title=element_blank())+
  coord_fixed(.7*max(ohc_intervals_df$year/3.5)/max(ohc_intervals_df$ohc_med))+
  ggtitle('Confidence and Credible Intervals for OHC')

if(!test_mode){
  save_image(pcredible,file=paste(figure_directory,'ohc_by_year.png',sep=''))
}else{
  plot(pcredible)
}

