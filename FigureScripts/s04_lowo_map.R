testmode=T
library(BayesianOHC)
library(tidyverse)
figure_directory='../Figures'

load('../ValidationData/lowodf_fullmodel.RData',verb=T)
grid_res=5

full_valid=lowodf_fullmodel%>%
  mutate(rmse=abs(validval-vhc_obs))

lowo_grid=full_valid%>%
  mutate(grid_lat=round((lat_degrees-grid_res/2)/grid_res)*grid_res+grid_res/2,
         grid_lon=round((lon_degrees-grid_res/2)/grid_res)*grid_res+grid_res/2)%>%
  group_by(grid_lat,grid_lon)%>%
  summarise(rmse=mean(rmse))%>%
  mutate(lat_degrees=grid_lat,lon_degrees=grid_lon)

breakswendpoints=c(seq(0,5,length=6),15)
varsym='RMSE'
mylabels= c(sapply(1:(length(breakswendpoints)-2),
                   function(ii) bquote(.(breakswendpoints[ii])~'\u2264'~
                                         .(varsym)~'<'~.(breakswendpoints[ii+1]))),
            bquote(.(breakswendpoints[length(breakswendpoints)-1])~'\u2264'~
                     .(varsym)))

p_lowo_full=plot_cut(lowo_grid,mybreaks=seq(0,5,length=6),varsym='s',
                     varname='rmse',legendname =bquote('GJ/m'^2),is_grid = T,mylabels=rev(mylabels))

if(!test_mode){
  save_image(p_lowo_full,
             file=paste(figure_directory,'validplot_lowo.png',sep=''))
}else{
  p_lowo_full
}

