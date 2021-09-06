testmode=T
library(BayesianOHC)
library(tidyverse)
figure_directory='../Figures'

hlines=c(-65,65)
vlines=c(-70,30,120)

boundlist=list()
boundlist[[1]]=rbind(c(-70,-30),
                     c(-70,-65),
                     c(30,-65),
                     c(30,-30),
                     c(-70,-30))%>%
  data.frame()%>%set_names('x','y')%>%mutate(group=1)
boundlist[[2]]=rbind(c(-100,20),
                     c(30,20),
                     c(30,-30),
                     c(-70,-30),
                     c(-70,0))%>%
  data.frame()%>%set_names('x','y')%>%mutate(group=2)
boundlist[[3]]=rbind(c(120,-65),
                     c(120,30))%>%
  data.frame()%>%set_names('x','y')%>%mutate(group=3)
boundlist[[4]]=rbind(c(120,20),
                     c(180,20))%>%
  data.frame()%>%set_names('x','y')%>%mutate(group=4)
boundlist[[5]]=rbind(c(-180,20),
                     c(-100,20))%>%
  data.frame()%>%set_names('x','y')%>%mutate(group=5)
boundaries=do.call(rbind,boundlist)

label_locs=rbind(data.frame(name='NA',fn='NA (North Atlantic)',x=-25,y=45),
                 data.frame(name='TA',fn='TA (Tropical Atlantic)',x=-25,y=0),
                 data.frame(name='SA',fn='SA (South Atlantic)',x=-25,y=-45),
                 data.frame(name='NWP',fn='NWP (North West Pacific)',x=-155,y=45),
                 data.frame(name='WTP',fn='WTP (West Tropical Atlantic)',x=-155,y=0),
                 data.frame(name='SP',fn='SP (South Pacific)',x=-155,y=-45),
                 data.frame(name='IN',fn='IN (Indian)',x=70,y=0),
                 data.frame(name='SI',fn='SI (South Indian)',x=70,y=-45),
                 data.frame(name='NEP',fn='NEP (North East Pacific)',x=165,y=45),
                 data.frame(name='NTP',fn='NTP (North Tropical Pacific)',x=165,y=0),
                 data.frame(name='SP',fn='SP (South Pacific)',x=165,y=-45))

p_definebasin=ggplot(mapping=aes(x=lon_degrees,y=lat_degrees))+
  scale_color_discrete()+
  geom_path(data=boundaries,aes(x=x,y=y,group=group),col='red')+
  geom_hline(aes(yintercept=c(-65,65,-30)),col='red')+
  borders("world",fill="black",colour="black")+
  coord_fixed(ratio = 1)+
  scale_y_continuous(expand = c(.01,.01),name='',
                     breaks=seq(-90,90,by=30),
                     labels=c('90S','60S','30S','0','30N','60N','90N'))+
  scale_x_continuous(expand = c(.01,.01),name='',
                     breaks=seq(-180,180,by=60),
                     labels=c(paste(-seq(-180,-60,by=60),'W'),'0',
                              paste(seq(60,180,by=60),'E')))+
  theme(legend.key.size=unit(.001,'cm'))+
  scale_color_discrete(name='Basin')+
  geom_point(data=label_locs,aes(x=x,y=y,col=fn), alpha = 0)+
  geom_text(data=label_locs,aes(x=x,y=y,label=name,col=fn),fontface='bold',cex=3)+
  guides(colour = guide_legend("Basin", override.aes = list(size = 1.5, alpha = 1)))

if(!testmode){
  save_image(p_definebasin,
             file=paste(figure_directory,'basin_definitions.png'))
}else{
  plot(p_definebasin)
}
