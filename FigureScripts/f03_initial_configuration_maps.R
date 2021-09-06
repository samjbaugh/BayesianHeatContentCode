testmode=T
figure_directory='../Figures'
load('../MCMC_Input/initial_parameters.RData',verb=T)

vis_grid=gen_masked_grid(1,1)
vis_grid=augment_data(vis_grid%>%mutate(years=2007,mu0=0,slope=0),initparams)

myplots=paramgrid_to_plots(vis_grid,configtype="Initialization",
                           which_plots = c('efflon','efflat','stnugget','stdev'))

if(!testmode){
  save_image(myplots$phi+ggtitle(''),file=paste(figure_directory,'globphi.png'))
  save_image(myplots$theta_lon+ggtitle(''),file=paste(figure_directory,'globlon.png'))
  save_image(myplots$theta_lat+ggtitle(''),file=paste(figure_directory,'globlat.png'))
  save_image(myplots$nugget+ggtitle(''),file=paste(figure_directory,'globnugget.png'))
}else{
  require('gridExtra')
  grid.arrange(myplots$phi,myplots$theta_lat,
               myplots$theta_lon,myplots$nugget,
               nrow=2,ncol=2)
}
