# 00_init_system_data.R: 
# Initializes system data for the BayesianOHC package.
# IN: RawData/mask.RData
# OUT: BayesOHC/sysdata.rda
# TEST: NA

load('../RawData/mask.RData',verb=T)
rg_mask=mygrid
default_linkfuns=list('phi'=exp,'nugget'=exp,
               'theta_lat'=exp,'theta_lon'=exp,
               'theta'=exp,'mu0'=identity,'slope'=identity)
default_invlinkfuns=list('phi'=log,'nugget'=log,
                  'theta_lat'=log,'theta_lon'=log,
                  'theta'=log,'mu0'=identity,'slope'=identity)
default_varname_list=c('phi','nugget','theta_lat',
                       'theta_lon','mu0','slope')

require('sf')
require('tmap')
require('maptools')
standard="+proj=longlat +datum=WGS84 +no_defs"
robin="+proj=robin +lon_0=-160 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs"

data(World)
wld.ll <- st_transform(World, standard)
wld.sp <- nowrapSpatialPolygons(as(wld.ll, "Spatial"), offset = 180-160)
wld.sf <- st_as_sf(wld.sp)
world.robin <- st_transform(wld.sf,robin)

usethis::use_data(rg_mask,default_linkfuns,default_invlinkfuns,
                  default_varname_list,standard,robin,world.robin,
                  internal = TRUE,overwrite = TRUE)
