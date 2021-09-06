testmode=T
library(BayesianOHC)
figure_directory='../Figures'
load('../MCMC_Output/Output_20k_run/sample_matrices.RData',verb=T)

#only use inits for knotpoints/hyperparameters which aren't sampled:
load('../MCMC_Input/initial_parameters.RData',verb=T)
hypers=initparams$hyperparam_list
mybasis=initparams$basis_fields[,1:4]

#coarse grid for evaluating parameter field means:
if(test_mode){
  gridlocs=gen_masked_grid(10,10)
}else{
  gridlocs=gen_masked_grid(2,2)
}
nbasis=dim(mybasis)[1]
ngrid=dim(gridlocs)[1]
lower_quantile_params=list('basis_fields'=mybasis,
                           'hyperparam_list'=hypers)
upper_quantile_params=list('basis_fields'=mybasis,
                           'hyperparam_list'=hypers)

for(varname in c('phi','nugget','theta_lat','theta_lon')){
  basis_param_mat=sample_matrices[[varname]]

  #convert from basis values to real values:
  myrange=hypers[[varname]][['range']]
  mu=hypers[[varname]][['mu']]
  phi=hypers[[varname]][['phi']]

  gridlocs=gridlocs%>%mutate(theta_lat=myrange,theta_lon=myrange)
  mybasis=mybasis%>%mutate(theta_lat=myrange,theta_lon=myrange)

  C_pred_basis=outer(1:ngrid,1:nbasis,cyl_cor_double,gridlocs,mybasis)
  C_basis_basis=outer(1:nbasis,1:nbasis,cyl_cor_double,mybasis,mybasis)
  chol_basis=chol(C_basis_basis+1e-7*diag(nbasis))

  varname_basis=paste('basis',varname,sep='_')
  realval_mat=default_linkfuns[[varname]](sqrt(phi)*
                     C_pred_basis%*%
                     backsolve(chol_basis,basis_param_mat)+mu)

  param_means=apply(realval_mat,2,mean)
  ind05=which(param_means==sort(param_means)[round(.05*length(param_means))])[1]
  ind95=which(param_means==sort(param_means)[round(.95*length(param_means))])[1]
  lower_quantile_params$basis_fields[[paste('basis',varname,sep='_')]]=
    basis_param_mat[,ind05]
  upper_quantile_params$basis_fields[[paste('basis',varname,sep='_')]]=
    basis_param_mat[,ind95]
}

#krig at higher resolution:
if(!test_mode){
  grid1x1=gen_masked_grid(1,1)
}else{
  grid1x1=gen_masked_grid(5,5)
}

lower_quantile_field=augment_data(grid1x1,lower_quantile_params,
                                  varname_list=c('theta_lat','theta_lon',
                                             'phi','nugget'))
upper_quantile_field=augment_data(grid1x1,upper_quantile_params,
                                  varname_list=c('theta_lat','theta_lon',
                                             'phi','nugget'))
myvarnames=c('stdev','stnugget','efflat','efflon')
plotvars=c('stdev','efflat','efflon')

plots05=paramgrid_to_plots(lower_quantile_field,configtype="Q05 Param Field",
                           which_plots=plotvars)
plots95=paramgrid_to_plots(upper_quantile_field,configtype="Q95 Param Field",
                           which_plots=plotvars)

for(varname in c('phi','theta_lat','theta_lon')){
  plots05[[varname]]=plots05[[varname]]+guides(fill=F,col=F)+ggtitle('')
  plots95[[varname]]=plots95[[varname]]+guides(fill=F,col=F)+ggtitle('')
}

if(!test_mode){
  for(varname in c('phi','theta_lat','theta_lon')){
    save_image(plots05[[varname]],file=
                 paste(figure_directory,varname,'_quant_low.png',sep=''))
    save_image(plots95[[varname]],file=
                 paste(figure_directory,varname,'_quant_high.png',sep=''))
  }
}else{
  require('gridExtra')
  grid.arrange(plots05$phi,plots95$phi,plots05$theta_lat,plots95$theta_lat,
         plots05$theta_lon,plots95$theta_lon,nrow=3,ncol=2)
}
