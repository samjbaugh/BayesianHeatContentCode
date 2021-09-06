#' Generates grid over the sphere using the mask saved in MCMC_Input/mask.RData
#' (for our purposes this mask is obtained from the Roemmich-Gilson climatology)
#'
#' @param latres The latitudinal resolution of the desired grid
#' @param lonres The longitudinal resolution of the desired grid
#' @return Returns a dataframe with columns "lat_rad", "lon_rad",
#' "lat_degrees", and "lon_degrees" corresponding to the gridpoint
#' locations
#' @export
gen_masked_grid<-function(latres=1,lonres=1)
{
  #load grid
  names(rg_mask)<-c("lon","lat","mask")

  #remove Mediterranean from mask
  in_med<-function(df)
  {
    return((df$lat>20 & df$lat<60) & (df$lon>2 & df$lon<45))
  }
  rg_mask$mask[in_med(rg_mask)]=0

  #coarsen grid
  newgrid=expand.grid(seq(min(rg_mask$lon),max(rg_mask$lon),by=lonres),
                      seq(min(rg_mask$lat),max(rg_mask$lat),by=latres))
  names(newgrid)<-c('lon','lat')
  numlons=length(unique(rg_mask$lon))
  masked_grid_degrees=newgrid[rg_mask$mask[1+round(newgrid$lon-min(rg_mask$lon))+
                                             numlons*round(newgrid$lat-min(rg_mask$lat))]==1,]
  masked_grid=data.frame('lat_rad'=deg2rad(masked_grid_degrees$lat),
                         'lon_rad'=deg2rad(masked_grid_degrees$lon),
                        'lat_degrees'=masked_grid_degrees$lat,
                        'lon_degrees'=masked_grid_degrees$lon)
  return(masked_grid)
}

