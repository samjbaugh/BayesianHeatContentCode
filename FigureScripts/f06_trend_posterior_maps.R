testmode=T
library(BayesianOHC)
figure_directory='../Figures'

load('../MCMC_Output/Output_20k_run/trend_resampled.RData')

pm_df=apply(all_field_samples,2,function(x) mean(x>0))

krig_grid=gen_masked_grid(1,1)
p05=plot_cut_diverging(krig_grid%>%mutate(slope=all_field_samples[low_quant_i,]),
                       seq(-.5,.5,length=6),
                       varsym=bquote(beta),varname='slope',
                       legendname=bquote((GJ/m^2)*'/year'),
                       mytitle='',base_map = mybasemap,is_grid=T)
a=round(all_int_means[low_quant_i]/10^12,3)
p05=p05+ggtitle(bquote('Q05 Integrated value:'~.(a)*'x'*10^21*'J/year'))

p95=plot_cut_diverging(krig_grid%>%mutate(slope=all_field_samples[up_quant_i,]),
                       seq(-.5,.5,length=6),
                       varsym=bquote(beta),varname='slope',
                       legendname=bquote((GJ/m^2)*'/year'),
                       mytitle='mytitle',base_map = mybasemap,is_grid=T)
a=round(intslope_samples[up_quant_i]/10^12,3)
p95=p95+ggtitle(bquote('Q95 Integrated value:'~.(a)*'x'*10^21*'J/year'))
p95

breakswendpoints=c(0,1,5,10,90,95,99,100)
varsym=bquote(beta)
mylabels= c(sapply(1:(length(breakswendpoints)-2),
                   function(ii) bquote(.(breakswendpoints[ii])*'% \u2264'~
                                         .(varsym)~'<'~.(breakswendpoints[ii+1])*'%')),
            bquote(.(breakswendpoints[length(breakswendpoints)-1])*'% \u2264'~
                     .(varsym)~'\u2264'~.(breakswendpoints[length(breakswendpoints)])*'%'))

breakswendpoints[1]=-1e-18
p_agreement=plot_cut_diverging(krig_grid%>%mutate(slope=pm_df*100),
                               breakswendpoints,varsym=bquote(beta),varname='slope',
                               legendname='',mytitle='',base_map = mybasemap,
                               is_grid=T,rounddigit = 1,
                               mylabels=mylabels)+
  ggtitle('Pointwise Posterior Probability of Positive Trend')


if(!test_mode){
  save_image(p05+ggtitle(''),
             file=paste(figure_directory,'slopeq05.png',sep=''))
  save_image(p_agreement+ggtitle(''),
             file=paste(figure_directory,'sign_agreement.png',sep=''))
  save_image(p95+ggtitle(''),
             file=paste(figure_directory,'slopeq95.png',sep=''))
}else{
  require('gridExtra')
  grid.arrange(p05,p_agreement,p95)
}
