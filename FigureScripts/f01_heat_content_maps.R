testmode=T
library(BayesianOHC)
figure_directory='../Figures'
load('../MCMC_Input/initial_parameters.RData',verb=T)

p_vhc_raw=p_vhc_raw_custom(augdata_init %>%
               subset(years==2016),mybreaks=seq(30,80,length=6),
             varname='vhc_obs',varsym='OHC',
             legendname=bquote("GJ/m"^2),mytitle='Vertical Column Heat Content (2016)',
             is_grid = F,mycex=.6) +
  guides(colour = guide_legend(bquote("GJ/m"^2), override.aes = list(size = 4)))+
  ggtitle('')

p_vhc_anom=p_vhc_anom_custom(augdata_init %>%
               subset(years==2016) %>% mutate(anom=vhc_obs-mu),
             mybreaks=seq(-3,3,length=6),
             varname='anom',varsym='Anom',
             legendname=bquote("GJ/m"^2),
             mytitle='Heat Content Anomalies (2016)',is_grid = F,mycex=.6) +
  guides(colour = guide_legend(bquote("GJ/m"^2), override.aes = list(size = 4)))+
  ggtitle('')

if(!testmode){
  save_image(p_vhc_raw,file=paste(figure_directory,
                                  'vhc_raw_2016.png'))
  save_image(p_vhc_anom,file=paste(figure_directory,
                                   'anom2016.png'))
}else{
  require('gridExtra')
  grid.arrange(p_vhc_raw,p_vhc_anom)
}
