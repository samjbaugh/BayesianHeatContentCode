#' Saves high-quality image and trims white space from result
#'
#' @param image Image to save
#' @param file Filename
#'
#' @import sf
#' @import tmap
#' @import maptools
#' @importFrom grDevices dev.off png
#' @importFrom methods as
#' @import ggplot2
#' @export
save_image<-function(image,file)
{
  grDevices::png(file, width=3200, height=3200, res=600)
  print(image)
  grDevices::dev.off()
  system(paste("convert -trim",file,file))
}

#' Plots map of observations "cut" by supplied breaks
#'
#' @param mydf Dataframe to plot
#' @param mybreaks Breakpoints to use
#' @param varsym Symbol of variablename for legend
#' @param varname Name of variable to plot
#' @param legendname Name of legend
#' @param mytitle Plot title
#' @param is_grid Is data located on grid?
#' @param rounddigit How many digits to round?
#' @param mycex Cex value for plotting if is_grid=F
#' @param mylabels Labels to use in legend
#' @export
plot_cut<-function(mydf,mybreaks,varsym,varname,
                   legendname,mytitle='',
                   is_grid=T,rounddigit=1,
                   mycex=1,mylabels=NULL)
{
  if(!is_grid){
    tt=mydf
    tt$x=tt$lon_degrees
    tt$y=tt$lat_degrees
    lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_point(cbind(tt$x[i],tt$y[i]))))

    breakswendpoints=mybreaks
    minval=min(mydf[[varname]])
    maxval=max(mydf[[varname]])
    if(min(breakswendpoints)>minval){
      breakswendpoints=c(minval-1e-4,breakswendpoints)
    }
    if(max(breakswendpoints)<maxval){
      breakswendpoints=c(breakswendpoints,maxval+1e-4)
    }
    mydf$varcut=cut(mydf[[varname]],c(breakswendpoints))
    print(breakswendpoints)
    breakswendpoints=sapply(breakswendpoints,function(x) round(x,rounddigit))

    pt0=st_sf(lll,varcut=mydf$varcut,crs=standard)
    pt2=st_transform(pt0,crs="+proj=robin +lon_0=-160 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs")

    if(is.null(mylabels)){
      mylabels=rev(sapply(1:(length(breakswendpoints)-1),function(ii)
        bquote(.(breakswendpoints[ii])~'\u2264'~.(varsym)~'<'~.(breakswendpoints[ii+1]))))
    }
    pcut=ggplot()+
      geom_sf(data = world.robin,col='black',fill='black')+
      geom_sf(data=pt2,aes(col=pt2$varcut,fill=pt2$varcut),cex=mycex)+
      # geom_point(data=knot_points,aes(x=lon_degrees,y=lat_degrees),col='black',size=.001,alpha=.5)+
      scale_color_viridis_d(name=legendname,drop=FALSE,direction=-1,
                            limits = rev(levels(mydf$varcut)),
                            labels=mylabels)+
      scale_fill_viridis_d(name=legendname,drop=FALSE,direction=-1,
                           limits = rev(levels(mydf$varcut)),
                           labels=mylabels)+
      ggtitle(mytitle)
  }else{
    tt=mydf
    tt$x=tt$lon_degrees#*6371*1000,
    tt$y=tt$lat_degrees#*6371*1000)

    dx=(min(diff(tt$lon_degrees)[diff(tt$lon_degrees)>0]))#*6371*1000
    dy=(min(diff(tt$lat_degrees)[diff(tt$lat_degrees)>0]))#*6371*1000

    lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_polygon(   list(rbind(c(tt$x[i],tt$y[i]),
                                                                        c(tt$x[i]+dx,tt$y[i]),
                                                                        c(tt$x[i]+dx,tt$y[i]+dy),
                                                                        c(tt$x[i],tt$y[i]+dy),
                                                                        c(tt$x[i],tt$y[i]))))))

    breakswendpoints=mybreaks
    minval=min(mydf[[varname]])
    maxval=max(mydf[[varname]])
    if(min(breakswendpoints)>minval){
      breakswendpoints=c(minval-1e-4,breakswendpoints)
    }
    if(max(breakswendpoints)<maxval){
      breakswendpoints=c(breakswendpoints,maxval+1e-4)
    }
    mydf$varcut=cut(mydf[[varname]],c(breakswendpoints))
    print(breakswendpoints)
    breakswendpoints=sapply(breakswendpoints,function(x) round(x,rounddigit))

    pt0=st_sf(lll,varcut=mydf$varcut,crs=standard)
    pt0.std=st_transform(pt0,crs=standard)
    pt0.sp <- nowrapSpatialPolygons(as_Spatial(pt0.std), offset = 180-160)
    pt0.sf <- st_as_sf(pt0.sp)
    pt0.sf$varcut=mydf$varcut
    pt0.robin <- st_transform(pt0.sf,robin)

    if(is.null(mylabels)){
      mylabels=rev( #c(bquote(.(varsym)~'<'~.(mybreaks[1])),
        sapply(1:(length(breakswendpoints)-1),function(ii) bquote(.(breakswendpoints[ii])~'\u2264'~.(varsym)~'<'~.(breakswendpoints[ii+1]))))
    }
    pcut=ggplot()+
      geom_sf(data = world.robin,col='grey',fill='grey')+
      geom_sf(data=pt0.robin,aes(col=pt0.robin$varcut,fill=pt0.robin$varcut))+
      # geom_point(data=knot_points,aes(x=lon_degrees,y=lat_degrees),col='black',size=.001,alpha=.5)+
      scale_color_viridis_d(name=legendname,drop=FALSE,direction=-1,
                            limits = rev(levels(mydf$varcut)),
                            labels=mylabels)+
      scale_fill_viridis_d(name=legendname,drop=FALSE,direction=-1,
                           limits = rev(levels(mydf$varcut)),
                           labels=mylabels)+
      coord_sf(expand=F,crs = robin)+ #ylim =  deg2rad(c(-80, 80))*6371*1000, xlim = deg2rad(c(-170, 170))*6371*1000, crs = robin)+
      ggtitle(mytitle)
  }
  pcut=pcut+
    theme(axis.text.x = element_text(size=7))
  return(pcut)
}

#' Plots map of observations "cut" by supplied breaks with diverging
#' colorbar
#'
#' @param mydf Dataframe to plot
#' @param mybreaks Breakpoints to use
#' @param varsym Symbol of variablename for legend
#' @param varname Name of variable to plot
#' @param legendname Name of legend
#' @param mytitle Plot title
#' @param is_grid Is data located on grid?
#' @param rounddigit How many digits to round?
#' @param mycex Cex value for plotting if is_grid=F
#' @param mylabels Labels to use in legend
#' @export
plot_cut_diverging<-function(mydf,mybreaks,varsym,varname,legendname,mytitle='',
                             is_grid=T,rounddigit=1,mycex=1,
                             mylabels=NULL)
{
  if(!is_grid){
    tt=mydf
    tt$x=tt$lon_rad*6371*1000
    tt$y=tt$lat_rad*6371*1000
    lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_point(cbind(tt$x[i],tt$y[i]))))

    breakswendpoints=mybreaks
    minval=min(mydf[[varname]])
    maxval=max(mydf[[varname]])
    if(min(breakswendpoints)>minval){
      breakswendpoints=c(minval-1e-4,breakswendpoints)
    }
    if(max(breakswendpoints)<maxval){
      breakswendpoints=c(breakswendpoints,maxval+1e-4)
    }
    mydf$varcut=cut(mydf[[varname]],c(breakswendpoints))
    print(breakswendpoints)
    breakswendpoints=sapply(breakswendpoints,function(x) round(x,rounddigit))

    pt0=st_sf(lll,varcut=mydf$varcut,crs=standard)
    pt2=st_transform(pt0,crs=robin)


    pcut=ggplot()+
      geom_sf(data = world.robin)+
      geom_sf(data=pt2,aes(col=pt2$varcut))+
      # geom_point(data=knot_points,aes(x=lon_degrees,y=lat_degrees),col='black',size=.001,alpha=.5)+
      scale_color_brewer(type='div',palette = "RdBu",name=legendname,drop=FALSE,direction=-1,
                         limits = (levels(mydf$varcut)),
                         labels=( #c(bquote(.(varsym)~'<'~.(mybreaks[1])),
                           sapply(1:(length(breakswendpoints)-1),function(ii) bquote(.(breakswendpoints[ii])~'\u2264'~.(varsym)~'<'~.(breakswendpoints[ii+1])))))+
      ggtitle(mytitle)+
      theme(plot.background  = element_rect(colour = "gray87", fill = "gray87"))
  }else{
    tt=mydf
    tt$x=tt$lon_degrees#*6371*1000,
    tt$y=tt$lat_degrees#*6371*1000)

    dx=(min(diff(tt$lon_degrees)[diff(tt$lon_degrees)>0]))#*6371*1000
    dy=(min(diff(tt$lat_degrees)[diff(tt$lat_degrees)>0]))#*6371*1000

    lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_polygon(   list(rbind(c(tt$x[i],tt$y[i]),
                                                                        c(tt$x[i]+dx,tt$y[i]),
                                                                        c(tt$x[i]+dx,tt$y[i]+dy),
                                                                        c(tt$x[i],tt$y[i]+dy),
                                                                        c(tt$x[i],tt$y[i]))))))

    breakswendpoints=mybreaks
    minval=min(mydf[[varname]])
    maxval=max(mydf[[varname]])
    if(min(breakswendpoints)>minval){
      breakswendpoints=c(minval-1e-4,breakswendpoints)
    }
    if(max(breakswendpoints)<maxval){
      breakswendpoints=c(breakswendpoints,maxval+1e-4)
    }
    mydf$varcut=cut(mydf[[varname]],c(breakswendpoints))
    print(breakswendpoints)

    pt0=st_sf(lll,varcut=mydf$varcut,crs=standard)
    pt0.std=st_transform(pt0,crs=standard)
    pt0.sp <- nowrapSpatialPolygons(as_Spatial(pt0.std), offset = 180-160)
    pt0.sf <- st_as_sf(pt0.sp)
    pt0.sf$varcut=mydf$varcut
    pt0.robin <- st_transform(pt0.sf,robin)

    breakswendpoints=sapply(breakswendpoints,function(x) round(x,rounddigit))

    ncols=length(levels(mydf$varcut))
    mycolvals=rev(RColorBrewer::brewer.pal(ncols,'RdBu'))
    mycolvals[ceiling(ncols/2)]='grey80'
    if(is.null(mylabels)){
      mylabels= c(sapply(1:(length(breakswendpoints)-2),
                         function(ii) bquote(.(breakswendpoints[ii])~'\u2264'~
                                               .(varsym)~'<'~.(breakswendpoints[ii+1]))),
                  bquote(.(breakswendpoints[length(breakswendpoints)-1])~'\u2264'~
                           .(varsym)~'\u2264'~.(breakswendpoints[length(breakswendpoints)])))
    }
    pcut=ggplot()+
      geom_sf(data = world.robin,col='black',fill='black')+
      geom_sf(data=pt0.robin,aes(col=pt0.robin$varcut,fill=pt0.robin$varcut))+
      # geom_point(data=knot_points,aes(x=lon_degrees,y=lat_degrees),col='black',size=.001,alpha=.5)+
      scale_color_manual(values=mycolvals,
                         name=legendname,drop=FALSE,
                         limits = (levels(mydf$varcut)),
                         labels=mylabels)+
      scale_fill_manual(values=mycolvals,
                        name=legendname,drop=FALSE,
                        limits = (levels(mydf$varcut)),
                        labels=mylabels)+#c(bquote(.(varsym)~'<'~.(mybreaks[1])),
      coord_sf(expand=F,crs = robin)+ #ylim =  deg2rad(c(-80, 80))*6371*1000, xlim = deg2rad(c(-170, 170))*6371*1000, crs = robin)+
      ggtitle(mytitle)


  }
  pcut=pcut+
    theme(axis.text.x = element_text(size=7))
  return(pcut)
}

#' Creates plots for each variable location on a grid
#'
#' @param plotgrid Gridded datframe containing values to plot
#' @param configtype Configuration type for labels
#' @param breaklist List of breakpoints for each variable
#' @param which_plots Which variables to plot
#' @export
paramgrid_to_plots<-function(plotgrid,configtype='',breaklist,
                             which_plots=c('stdev','nugget','efflat',
                                           'efflon','mu','slope')){
  is_grid=T
  if(missing(breaklist)){
    breaklist=list()
    breaklist[['efflat']]=exp(seq(0,log(40),length=6))
    breaklist[['nugget']]=seq(0,.5,length=6)
    breaklist[['stdev']]=exp(seq(sqrt(1),log(sqrt(90)),length=6))
    breaklist[['efflon']]=exp(seq(0,log(40),length=6))
    breaklist[['mu']]=seq(0,80,length=10)
    breaklist[['slope']]= seq(-.5,.5,length=10)
  }

  plotgrid$efflat=convert_theta_lat_to_effective_range_deg(plotgrid$theta_lat)
  plotgrid$efflon=convert_theta_lon_to_effective_range_deg(plotgrid$theta_lon)
  plotgrid$stdev=sqrt(plotgrid$phi)
  plotgrid$stnugget=sqrt(plotgrid$nugget*plotgrid$phi)

  varname='stdev'
  if('stdev'%in%which_plots)
  {
    mybreaks=exp(seq(log(sqrt(1)),log(sqrt(90)),length=6))
    varsym=bquote(sqrt(varphi))
    mytitle=paste('Marginal Standard Deviation (',configtype,')',sep='') #bquote(sqrt(varphi) ~ .(configtype))
    pstdev=plot_cut(plotgrid,mybreaks,varsym=varsym,
                    varname=varname,legendname=bquote('GJ/m'^2),mytitle=mytitle,
                    is_grid=is_grid)
    pstdev
  }else{
    pstdev=c()
  }


  #nugget
  varname='stnugget'
  if('stnugget'%in%which_plots)
  {
    plotgrid$efflat=convert_theta_lat_to_effective_range_deg(plotgrid$theta_lat)
    plotgrid$efflon=convert_theta_lon_to_effective_range_deg(plotgrid$theta_lon)
    plotgrid$stdev=sqrt(plotgrid$phi)
    plotgrid$stnugget=sqrt(plotgrid$nugget*plotgrid$phi)
    varsym=bquote(sigma)
    mytitle=bquote(sigma ~ .(configtype))
    # mybreaks=seq(exp(-10),exp(1),length=6)
    mybreaks=exp(seq(log(sqrt(1)),log(sqrt(5)),length=6))
    mybreaks=seq(0,2,length=6)
    mybreaks=exp(seq(-1,log(2),length=6))
    mytitle='' #paste('Standard Deviation of Nugget (',configtype,')',sep='') #bquote(sqrt(varphi) ~ .(configtype))
    pnugget=plot_cut(plotgrid,mybreaks,varsym=varsym,
                     varname=varname,legendname=bquote('GJ/m'^2),mytitle=mytitle,
                     is_grid=is_grid)
  }else{
    pnugget=c()
  }

  #effective lat range
  varname='efflat'
  if('efflat'%in%which_plots)
  {
    mybreaks=exp(seq(1,log(40),length=6))
    varsym=bquote(gamma[lat])
    mytitle=paste('Effecitve Latitudinal Range (',configtype,')',sep='') #bquote(sqrt(varphi) ~ .(configtype))
    pefflat=plot_cut_theta_lat_custom(plotgrid,mybreaks,varsym=varsym,
                                      varname=varname,legendname='Degrees Latitude',mytitle=mytitle,
                                      is_grid=is_grid)
  }else{
    pefflat=c()
  }

  #effective lon range
  varname='efflon'
  if('efflon'%in%which_plots)
  {
    mybreaks=exp(seq(1,log(40),length=6))
    varsym=bquote(gamma[lon])
    mytitle=paste('Effecitve Longitudinal Range (',configtype,')',sep='') #bquote(sqrt(varphi) ~ .(configtype))
    pefflon=plot_cut(plotgrid,mybreaks,varsym=varsym,
                     varname=varname,legendname='Degrees Longitude',mytitle=mytitle,
                     is_grid=is_grid)
  }else{
    pefflon=c()
  }

  varname='mu'
  if('mu'%in%which_plots){
    mybreaks=seq(30,80,length=6)
    varsym=bquote(mu[2007])
    mytitle=bquote(mu ~ .(configtype))
    mytitle=paste('Mean Field for 2007 (',configtype,')',sep='') #bquote(sqrt(varphi) ~ .(configtype))
    pmu=plot_cut(plotgrid,mybreaks,varsym=varsym,
                 varname=varname,legendname=bquote('GJ/m'^2),mytitle=mytitle,
                 is_grid=is_grid)
  }else{
    pmu=c()
  }
  varname='slope'
  if('slope'%in%which_plots)
  {
    mybreaks=seq(-.5,.5,length=6)
    varname='slope'
    varsym=bquote(beta)
    mytitle=bquote(beta ~ .(configtype))
    mytitle=paste('Yearly Trend (',configtype,')') #bquote(sqrt(varphi) ~ .(configtype))
    pslope=plot_cut_diverging(plotgrid,mybreaks,varsym=varsym,
                              varname=varname,legendname=bquote('GJ/(m'^2*'year'),
                              mytitle=mytitle,is_grid=is_grid)
  }else{
    pslope=c()
  }

  return(list('phi'=pstdev,'theta_lon'=pefflon,'theta_lat'=pefflat,
              'nugget'=pnugget,'mu'=pmu,'slope'=pslope))
}

p_vhc_raw_custom<-function(mydf,mybreaks,varsym,varname,legendname,mytitle='',
                           is_grid=T,rounddigit=1,
                           mycex=1,mylabels=NULL,div=F)
{
  tt=mydf
  tt$x=tt$lon_degrees
  tt$y=tt$lat_degrees
  lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_point(cbind(tt$x[i],tt$y[i]))))

  breakswendpoints=mybreaks
  minval=min(mydf[[varname]])
  maxval=max(mydf[[varname]])
  if(min(breakswendpoints)>minval){
    breakswendpoints=c(minval-1e-4,breakswendpoints)
  }
  if(max(breakswendpoints)<maxval){
    breakswendpoints=c(breakswendpoints,maxval+1e-4)
  }
  mydf$varcut=cut(mydf[[varname]],c(breakswendpoints))
  print(breakswendpoints)
  breakswendpoints=sapply(breakswendpoints,function(x) round(x,rounddigit))

  pt0=st_sf(lll,varcut=mydf$varcut,crs=standard)
  pt2=st_transform(pt0,crs="+proj=robin +lon_0=-160 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs")

  if(is.null(mylabels)){
    mylabels=rev(sapply(1:(length(breakswendpoints)-1),function(ii)
      bquote(.(breakswendpoints[ii])~'\u2264'~.(varsym)~'<'~.(breakswendpoints[ii+1]))))
  }


  tt=gen_masked_grid(2,2)
  tt$x=tt$lon_degrees#*6371*1000,
  tt$y=tt$lat_degrees#*6371*1000)

  dx=(min(diff(tt$lon_degrees)[diff(tt$lon_degrees)>0]))#*6371*1000
  dy=(min(diff(tt$lat_degrees)[diff(tt$lat_degrees)>0]))#*6371*1000

  lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_polygon(   list(rbind(c(tt$x[i],tt$y[i]),
                                                                      c(tt$x[i]+dx,tt$y[i]),
                                                                      c(tt$x[i]+dx,tt$y[i]+dy),
                                                                      c(tt$x[i],tt$y[i]+dy),
                                                                      c(tt$x[i],tt$y[i]))))))

  pt3=st_sf(lll,varcut=mydf$varcut,crs=standard)
  pt3.std=st_transform(pt3,crs=standard)
  pt3.sp <- nowrapSpatialPolygons(as_Spatial(pt3.std),
                                  offset = 180-160)
  pt3.sf <- st_as_sf(pt3.sp)
  pt3.robin <- st_transform(pt3.sf,robin)

  pcut=ggplot()+
    geom_sf(data = world.robin,col='black',fill='black')+
    geom_sf(data=pt3.robin,fill='darkgrey',col='darkgrey',cex=mycex)+
    geom_sf(data=pt2,aes(col=pt2$varcut),cex=mycex)
  if(div){
    pcut=pcut+scale_color_brewer(type='div',palette = "RdBu",name=legendname,drop=FALSE,direction=-1,
                                 limits = rev(levels(mydf$varcut)),
                                 labels=mylabels)
  }else{
    pcut=pcut+scale_color_viridis_d(name=legendname,drop=FALSE,direction=-1,
                                    limits = rev(levels(mydf$varcut)),
                                    labels=mylabels)+
      ggtitle(mytitle)
  }
  pcut=pcut+
    theme(axis.text.x = element_text(size=7))+
    coord_sf(expand=F,crs = robin)
  return(pcut)
}


p_vhc_anom_custom<-function(mydf,mybreaks,varsym,varname,legendname,mytitle='',
                            is_grid=T,rounddigit=1,
                            mycex=1,mylabels=NULL,div=F)
{
  tt=mydf
  tt$x=tt$lon_degrees
  tt$y=tt$lat_degrees
  lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_point(cbind(tt$x[i],tt$y[i]))))

  breakswendpoints=mybreaks
  minval=min(mydf[[varname]])
  maxval=max(mydf[[varname]])
  if(min(breakswendpoints)>minval){
    breakswendpoints=c(minval-1e-4,breakswendpoints)
  }
  if(max(breakswendpoints)<maxval){
    breakswendpoints=c(breakswendpoints,maxval+1e-4)
  }
  mydf$varcut=cut(mydf[[varname]],c(breakswendpoints))
  print(breakswendpoints)
  breakswendpoints=sapply(breakswendpoints,function(x) round(x,rounddigit))

  pt0=st_sf(lll,varcut=mydf$varcut,crs=standard)
  pt2=st_transform(pt0,crs="+proj=robin +lon_0=-160 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs")

  if(is.null(mylabels)){
    mylabels=rev(sapply(1:(length(breakswendpoints)-1),function(ii)
      bquote(.(breakswendpoints[ii])~'\u2264'~.(varsym)~'<'~.(breakswendpoints[ii+1]))))
  }
  ncols=length(levels(mydf$varcut))
  mycolvals=(RColorBrewer::brewer.pal(ncols,'RdBu'))
  # mycolvals=c('#762a83','#af8dc3','#e7d4e8','#f7f7f7','#d9f0d3','#7fbf7b','#1b7837')
  mycolvals[ceiling(ncols/2)]='grey80'
  pcut=ggplot()+
    geom_sf(data = world.robin,col='black',fill='black')+
    geom_sf(data=pt2,aes(col=pt2$varcut),cex=mycex)+
    scale_color_manual(values=(mycolvals),
                       name=legendname,drop=FALSE,
                       limits = rev(levels(mydf$varcut)),
                       labels=(mylabels))+
    theme(axis.text.x = element_text(size=7))+
    coord_sf(expand=F,crs = robin)
  return(pcut)
}

plot_cut_theta_lat_custom<-function(mydf,mybreaks,varsym,varname,legendname,mytitle='',
                                    is_grid=T,rounddigit=1,
                                    mycex=1,mylabels=NULL)
{
  tt=mydf
  tt$x=tt$lon_degrees#*6371*1000,
  tt$y=tt$lat_degrees#*6371*1000)

  dx=(min(diff(tt$lon_degrees)[diff(tt$lon_degrees)>0]))#*6371*1000
  dy=(min(diff(tt$lat_degrees)[diff(tt$lat_degrees)>0]))#*6371*1000

  lll=sapply(1:dim(tt)[1],function(i) st_sfc(st_polygon(   list(rbind(c(tt$x[i],tt$y[i]),
                                                                      c(tt$x[i]+dx,tt$y[i]),
                                                                      c(tt$x[i]+dx,tt$y[i]+dy),
                                                                      c(tt$x[i],tt$y[i]+dy),
                                                                      c(tt$x[i],tt$y[i]))))))

  breakswendpoints=mybreaks
  minval=min(mydf[[varname]])
  maxval=max(mydf[[varname]])
  if(min(breakswendpoints)>minval){
    breakswendpoints=c(minval-1e-4,breakswendpoints)
  }
  if(max(breakswendpoints)<maxval){
    breakswendpoints=c(breakswendpoints,maxval+1e-4)
  }
  mydf$varcut=cut(mydf[[varname]],c(breakswendpoints))
  levels(mydf$varcut)=c(levels(mydf$varcut),'(40,Inf]')
  print(breakswendpoints)
  breakswendpoints=sapply(breakswendpoints,function(x) round(x,rounddigit))

  pt0=st_sf(lll,varcut=mydf$varcut,crs=standard)
  pt0.std=st_transform(pt0,crs=standard)
  pt0.sp <- nowrapSpatialPolygons(as(pt0.std, "Spatial"), offset = 180-160)
  pt0.sf <- st_as_sf(pt0.sp)
  pt0.sf$varcut=mydf$varcut
  pt0.robin <- st_transform(pt0.sf,robin)

  if(is.null(mylabels)){
    mylabels=rev( #c(bquote(.(varsym)~'<'~.(mybreaks[1])),
      sapply(1:(length(breakswendpoints)-1),function(ii) bquote(.(breakswendpoints[ii])~'\u2264'~.(varsym)~'<'~.(breakswendpoints[ii+1]))))
  }
  mylabels=c(bquote(.('40')~'\u2264'~.(varsym)),mylabels)
  pcut=ggplot()+
    geom_sf(data = world.robin,col='grey',fill='grey')+
    geom_sf(data=pt0.robin,aes(col=pt0.robin$varcut,fill=pt0.robin$varcut))+
    # geom_point(data=knot_points,aes(x=lon_degrees,y=lat_degrees),col='black',size=.001,alpha=.5)+
    scale_color_viridis_d(name=legendname,drop=FALSE,direction=-1,
                          limits = rev(levels(mydf$varcut)),
                          labels=mylabels)+
    scale_fill_viridis_d(name=legendname,drop=FALSE,direction=-1,
                         limits = rev(levels(mydf$varcut)),
                         labels=mylabels)+
    coord_sf(expand=F,crs = robin)+ #ylim =  deg2rad(c(-80, 80))*6371*1000, xlim = deg2rad(c(-170, 170))*6371*1000, crs = robin)+
    ggtitle(mytitle)
  pcut=pcut+
    theme(axis.text.x = element_text(size=7))
  return(pcut)
}

