#'Computes euclidean correlations using kernel convolutions
#' @description Computes the exact latitudinal correlation kernel convolution between
#' locations loc1 and loc2 with latitudinal ranges range1 and range2 respectively
#' @param loc1 The first location on the longitudinal circle
#' @param loc2 The second location on the longitudinal circle
#' @param range1 The longitudinal range for the first location
#' @param range2 The longitudinal range for the second location
#' @return Returns the exact value of of the convolution in the Euclidean (latitudinal)
#' dimension
#' @export
euclidean_correlation_convolution<-function(loc1,loc2,range1,range2)
{
  q_term=(loc1-loc2)^2/((range1+range2)/2)
  param_term=sqrt(2)*(range1^(1/4))*(range2^(1/4))/sqrt(range1+range2)
  retval=param_term*exp(-q_term)
  return(retval)
}

#' Computes cylindrical correlations using kernel convolutions
#' @description Computes the exact longitudinal correlation kernel convolution between
#' locations loc1 and loc2 with longitudinal ranges range1 and range2
#' respectively. Since this is over a circular domain Gaussian error function
#' calls are required.
#' @param loc1 The first location on the longitudinal circle
#' @param loc2 The second location on the longitudinal circle
#' @param range1 The longitudinal range for the first location
#' @param range2 The longitudinal range for the second location
#' @return Returns the value of the exact cylndrical correlation; see supplementary
#' material for details
#' @export
cylindrical_correlation_exact<-function(loc1,loc2,range1,range2){
  yprime=mygc_dist(loc1,loc2)
  term1=gaussian_integral(0,yprime,range1,range2,yprime-pi,pi)
  term2=gaussian_integral(2*pi,yprime,range1,range2,pi,yprime+pi)
  scaling_term=2/(sqrt(pi)*sqrt(sqrt(range1*range2)))
  retval=(term1+term2)*scaling_term
  return(retval)
}

#' Internal function for computing
#' exact kernel convolutions over the cylindrical dimension.
#'
#' @param loc1 First location
#' @param loc2 Second location
#' @param lower Lower bound
#' @param upper Upper bound
#' @param range1 Range parameter for first location
#' @param range2 Range parameter for second location
gaussian_integral<-function(loc1,loc2,range1,range2,upper,lower){
  d=(loc1-loc2)^2
  sdterm=(1/2)*(sqrt(range1*range2)/sqrt(range1+range2))
  tt=pnorm(upper,mean=(range2*loc1+range1*loc2)/(range1+range2),sd=sdterm)-
    pnorm(lower,mean=(range2*loc1+range1*loc2)/(range1+range2),sd=sdterm)
  return(-sqrt(2*pi)*sdterm*tt*exp(-d/((range1+range2)/2)))
}


#' Cylindrical correlations using Gaussian approximation
#' @description Computes the approximate longitudinal correlation kernel convolution between
#' locations loc1 and loc2 with longitudinal ranges range1 and range2
#' respectively.
#'
#' @param loc1 The first location on the longitudinal circle
#' @param loc2 The second location on the longitudinal circle
#' @param range1 The longitudinal range for the first location
#' @param range2 The longitudinal range for the second location
#' Returns the Gaussian approximation to the full convolution, which is accurate
#' if the effective range is less than about 100 degrees; see supplementary
#' material for details.
#' @export
cylindrical_correlation_gaussian<-function(loc1,loc2,range1,range2)
{
  myd=mygc_dist(loc1,loc2)
  q_term=myd^2/((range1+range2)/2)
  param_term=sqrt(2)*(range1^(1/4))*(range2^(1/4))/sqrt(range1+range2)
  retval=param_term*exp(-q_term)
  return(retval)
}

#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param myinput First dataframe
#'@export
cyl_cor_single=function(ii,jj,myinput) {
  (euclidean_correlation_convolution(myinput$lat_rad[ii],myinput$lat_rad[jj],
                        myinput$theta_lat[ii],myinput$theta_lat[jj])*
    cylindrical_correlation_gaussian(myinput$lon_rad[ii],myinput$lon_rad[jj],
                            myinput$theta_lon[ii],
                            myinput$theta_lon[jj]) +
    (ii==jj)*myinput$nugget[ii])
}

#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param myinput First dataframe
#'@export
cyl_cor_convolution_single=function(ii,jj,myinput) {
  (euclidean_correlation_convolution(myinput$lat_rad[ii],myinput$lat_rad[jj],
                                     myinput$theta_lat[ii],myinput$theta_lat[jj])*
     cylindrical_correlation_exact(myinput$lon_rad[ii],myinput$lon_rad[jj],
                                    myinput$theta_lon[ii],
                                    myinput$theta_lon[jj]) +
     (ii==jj)*myinput$nugget[ii])
}


#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param input1 First dataframe
#'@param input2 Second dataframe
#'@export
cyl_cor_double=function(ii,jj,input1,input2) {
  euclidean_correlation_convolution(input1$lat_rad[ii],input2$lat_rad[jj],
                        input1$theta_lat[ii],input2$theta_lat[jj])*
    cylindrical_correlation_gaussian(input1$lon_rad[ii],input2$lon_rad[jj],
                            input1$theta_lon[ii],input2$theta_lon[jj])
}

#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param input1 First dataframe
#'@param input2 Second dataframe
#'@export
cyl_cor_exact_double=function(ii,jj,input1,input2) {
  euclidean_correlation_convolution(input1$lat_rad[ii],input2$lat_rad[jj],
                                    input1$theta_lat[ii],input2$theta_lat[jj])*
    c(input1$lon_rad[ii],input2$lon_rad[jj],
                                     input1$theta_lon[ii],input2$theta_lon[jj])
}

#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param myinput First dataframe
#'@export
cyl_cor_exact_single=function(ii,jj,myinput) {
  (euclidean_correlation_convolution(myinput$lat_rad[ii],myinput$lat_rad[jj],
                                     myinput$theta_lat[ii],myinput$theta_lat[jj])*
     cylindrical_correlation_exact(myinput$lon_rad[ii],myinput$lon_rad[jj],
                                      myinput$theta_lon[ii],
                                      myinput$theta_lon[jj]) +
     (ii==jj)*myinput$nugget[ii])
}


#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param myinput First dataframe
#'@export
chord_cor_single=function(ii,jj,myinput) {
  range1=myinput$theta[ii]
  range2=myinput$theta[jj]

  distsq=(myinput$coord1[ii]-myinput$coord1[jj])^2+
    (myinput$coord2[ii]-myinput$coord2[jj])^2+
    (myinput$coord3[ii]-myinput$coord3[jj])^2

  q_term=distsq/((range1+range2)/2)
  param_term=(sqrt(2)*(range1^(1/4))*(range2^(1/4))/sqrt(range1+range2))^3
  retval=param_term*exp(-q_term)
  return(retval)
}

#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param input1 First location dataframe (needs coord1,coord2,coord3)
#'@param input2 Second location dataframe (needs coord1,coord2,coord3)
#'@export
chord_cor_double=function(ii,jj,input1,input2) {
  range1=input1$theta[ii]
  range2=input2$theta[jj]

  distsq=(input1$coord1[ii]-input2$coord1[jj])^2+
    (input1$coord2[ii]-input2$coord2[jj])^2+
    (input1$coord3[ii]-input2$coord3[jj])^2

  q_term=distsq/((range1+range2)/2)
  param_term=(sqrt(2)*(range1^(1/4))*(range2^(1/4))/sqrt(range1+range2))^3
  retval=param_term*exp(-q_term)
  return(retval)
}


#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param myinput Location dataframe (needs coord1,coord2,coord3)
#'@export
chord_cor_exp_single=function(ii,jj,myinput) {
  range1=myinput$theta[ii]
  range2=myinput$theta[jj]

  distsq=(myinput$coord1[ii]-myinput$coord1[jj])^2+
    (myinput$coord2[ii]-myinput$coord2[jj])^2+
    (myinput$coord3[ii]-myinput$coord3[jj])^2

  q_term=sqrt(distsq)/((range1+range2)/2)
  param_term=(sqrt(2)*(range1^(1/4))*(range2^(1/4))/sqrt(range1+range2))^3
  retval=param_term*exp(-q_term)
  return(retval)
}

#' Helper function for computing correlations
#'@param ii First index
#'@param jj Second index
#'@param input1 First location dataframe (needs coord1,coord2,coord3)
#'@param input2 Second location dataframe (needs coord1,coord2,coord3)
#'@export
chord_cor_exp_double=function(ii,jj,input1,input2) {
  range1=input1$theta[ii]
  range2=input2$theta[jj]

  distsq=(input1$coord1[ii]-input2$coord1[jj])^2+
    (input1$coord2[ii]-input2$coord2[jj])^2+
    (input1$coord3[ii]-input2$coord3[jj])^2

  q_term=sqrt(distsq)/((range1+range2)/2)
  param_term=(sqrt(2)*(range1^(1/4))*(range2^(1/4))/sqrt(range1+range2))^3
  retval=param_term*exp(-q_term)
  return(retval)
}
