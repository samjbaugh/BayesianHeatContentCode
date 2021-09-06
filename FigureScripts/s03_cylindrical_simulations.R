testmode=T
library(BayesianOHC)
library(tidyverse)
figure_directory='../Figures'

load('../ValidationData/conv_test_simulations.RData',verb=T)
load('../MCMC_Output/Output_20k_run/map_configuration.RData',verb=T)

q99=augment_data(gen_masked_grid(2,2),
                      map_sample$stored_params,varname_list='theta_lon')%>%
  .$theta_lon%>%convert_theta_lat_to_effective_range_deg%>%quantile(.99)

pcylind_sim=ggplot()+
  geom_line(data=stationary_simulations,aes(x=true_effrange,
                                            y=(approx-conv)/conv,col='Stationary'))+
  geom_line(data=nonstationary_simulations,aes(x=true_effrange,
                                               y=(approx-conv)/conv,col='Nonstationary'))+
  scale_x_continuous(name='True Effective Range',
                     breaks=seq(0,125,by=15),
                     labels=sapply(seq(0,125,by=15),function(x) paste(x,'°',sep='')))+
  scale_y_continuous(name='Fractional Error',
                     breaks=seq(-.1,0,length=6),
                     labels=sapply(seq(-.1,0,length=6),
                                   function(x) paste(x*100,'%',sep='')))+
  # geom_line(aes(x=q99,y=seq(-.1,.001,length=100),col='99% Range Quantile'))
  geom_vline(aes(xintercept=q99,col='99% Range Quantile'),show.legend = F)+
  scale_color_discrete(breaks=c('Stationary','Nonstationary','99% Range Quantile'),
                       labels=c('Stationary  ','Nonstationary  ',
                                bquote('99%'~theta[lat]~'init')))+
  theme(legend.title=element_blank())+
  coord_fixed(500)+
  theme(title=element_text(size=8),
        axis.text=element_text(size=6),
        axis.title=element_text(size=8),
        legend.text = element_text(size=6),
        legend.title=element_blank(),
        legend.margin =margin(c(-1,0,0,0)))
pcylind_sim
stationary_simulations%>%mutate(fracerror=(approx-conv)/conv)%>%
  dplyr::filter(abs(true_effrange-q99)==min(abs(true_effrange-q99)))%>%.$fracerror
nonstationary_simulations%>%mutate(fracerror=(approx-conv)/conv)%>%
  dplyr::filter(abs(true_effrange-q99)==min(abs(true_effrange-q99)))%>%.$fracerror

if(!test_mode){
  save_image(pcylind_sim,
             file=paste(figure_directory,'cylind_sim.png',sep=''))
}

