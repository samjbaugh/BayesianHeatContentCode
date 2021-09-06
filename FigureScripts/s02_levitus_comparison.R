testmode=T
library(BayesianOHC)
library(tidyverse)
figure_directory='../Figures'

load('../ValidationData/levitus_ohc.RData',verb=T)
load('../MCMC_Output/Output_20k_run/posterior_intervals.RData',verb=T)
load('../MCMC_Output/Output_20k_run/trend_resamp_intonly.RData',verb=T)

zalpha=qnorm(.975)

adj=.1
ebarwidth=.2
plev_comparison=ggplot(data=ohc_intervals_df)+
  geom_point(aes(x=year-adj,y=ohc_med,col='BayesianOHC'),cex=.7)+
  geom_point(data=levitus_ohc,aes(x=year+adj,y=ohc,col='Levitus et al.'),cex=.7)+
  geom_errorbar(aes(x=year-adj,ymin=lower_cred,ymax=upper_cred,col='BayesianOHC'),
                stat="identity",alpha=.9,width=ebarwidth)+
  geom_errorbar(data=levitus_ohc,aes(x=year+adj,
                                     ymin=ohc-zalpha*ohc_sd,
                                     ymax=ohc+zalpha*ohc_sd,col='Levitus et al.'),
                stat="identity",alpha=.9,width=ebarwidth)+
  scale_y_continuous(name=bquote('OHC (10'^25*' Joules)'),
                     breaks=seq(1513e22,1531e22,by=4e22)/1e9,
                     labels=seq(1.513,1.531,by=.004))+
  scale_x_continuous(breaks=seq(2007,2016,by=1),expand=c(.04,.04),name='Year')+
  scale_color_manual(values=c('red','blue'))+
  theme(legend.title=element_blank())+
  coord_fixed(.7*max(ohc_intervals_df$year/4.5)/max(ohc_intervals_df$ohc_med))+
  # ggtitle('OHC Uncertainty Intervals')+
  # guides(shape = guide_legend(override.aes = list(size = .3)),
         # color = guide_legend(override.aes = list(size = .3))) +
  theme(title=element_text(size=8),
        axis.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.text = element_text(size=6),
        legend.title=element_blank(),
        legend.margin =margin(c(-1,-6,0,-6)))

plev_comparison

if(!test_mode){
  save_image(plev_comparison,
             file=paste(figure_directory,'levitus_comparison.png',sep=''))
}else{
  plot(plev_comparison)
}


