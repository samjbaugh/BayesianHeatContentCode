#' Computes CRPS values
#' @param x Validation dataframe
#' @importFrom verification crps
#' @export
get_crps<-function(x){verification::crps(x$vhc_obs,cbind(x$validval,x$validsd))$CRPS}

#' Computes RMSE values
#' @param x Validation dataframe
#' @export
get_rmse<-function(x){mean(abs(x$validval-x$vhc_obs))}

#' Adds 3D coordinates to dataframe
#'
#' @param mydf Input dataframe
#' @param R Radius of the earth; defaults to unit circle but uses R=6361 if
#' "earth" is used as the argument
#' @return Returns mydf with coord1, coord2, and coord3 added
#' @export
add_3d_coord<-function(mydf,R=1){
  if(R=='earth'){
    R=6371
  }
  cbind(mydf,lat_lon_to_chordal(mydf$lat_rad,mydf$lon_rad,R=R))
}

#' Generates delta values on a grid for computing numerical integrals
#'
#' @param grid_locs Grid locations
#' @param res Resolution of the grid; if not included will be inferred
#' @export
gen_deltas_from_grid<-function(grid_locs,res=NULL){
  if(is.null(res)){
    res=diff(grid_locs$lon_rad)[1]
  }
  circ=40075*10^3 #circumference in meters
  latres_dense_m=res*circ/(2*pi)
  lonres_dense_m=res*cos(grid_locs$lat_rad)*circ/(2*pi)
  scalarvec=latres_dense_m*lonres_dense_m
  return(scalarvec)
}

#' Converts theta lon to effective range in degrees
#'
#' @param thetalon Theta lon
#' @export
convert_theta_lon_to_effective_range_deg<-function(thetalon)
{
  delta_rad_lon=sqrt(-log(0.05)*(thetalon))
  r=6371
  return(rad2deg(delta_rad_lon))
}

#' Converts theta lat to effective range in degrees
#'
#' @param thetalat Theta lat
#' @export
convert_theta_lat_to_effective_range_deg<-function(thetalat)
{
  delta_rad_lat=sqrt(-log(0.05)*(thetalat))
  r=6371
  return(rad2deg(delta_rad_lat))
}

#' Compute great circle distance in radians
#'
#' @param l1 Argument 1
#' @param l2 Argument 2
#' @export
mygc_dist<-function(l1,l2)
{
  return(pmin(abs(l1-l2),2*base::pi-abs(l1-l2)))
}

#' Compute great circle distance in degrees
#'
#' @param l1 Argument 1
#' @param l2 Argument 2
#' @export
mygc_dist_degrees<-function(l1,l2)
{
  return(pmin(abs(l1-l2),360-abs(l1-l2)))
}


lat_lon_to_chordal<-function(mylat,mylon,R=NULL){
  if(is.null(R)){
    R=6371 # efault to Earth mean radius [km]
  }
  x = R*cos(mylon)*sin(mylat)
  y = R*sin(mylon)*sin(mylat)
  z = R*cos(mylat)
  return(data.frame('coord1'=x,'coord2'=y,'coord3'=z))
}

rad2deg<-function(rad)
{
  return(rad*180/base::pi)
}

deg2rad<-function(deg)
{
  return(deg*base::pi/180)
}

chol_solve<-function(mychol,myvec)
{
  return(backsolve(mychol,forwardsolve(t(mychol),myvec)))
}

myeuc_dist<-function(l1,l2)
{
  return(abs(l1-l2))
}

window_dist<-function(x1,x2)
{
  gg=function(ii,jj) pmax(abs(x1$lat_rad[ii]-x2$lat_rad[jj]),
                          mygc_dist(x1$lon_rad[ii],x2$lon_rad[jj]))
  n1=dim(x1)[1]
  n2=dim(x2)[1]
  return(outer(1:n1,1:n2,gg))
}

gen_pdist_mat_cylindrical<-function(x1,x2)
{
  gg=function(ii,jj) sqrt((x1$lat_rad[ii]-x2$lat_rad[jj])^2+
                            mygc_dist(x1$lon_rad[ii],x2$lon_rad[jj])^2)
  n1=dim(x1)[1]
  n2=dim(x2)[1]
  return(outer(1:n1,1:n2,gg))
}





