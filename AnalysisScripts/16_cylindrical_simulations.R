# 16_cylindrical_simulations:
# Runs the simulation study for comparing the accuracy of the exact
# cylindrical convolutions versus the Gaussian approximation; see Section S3
# in the supplementary materials for details. Most of the work is done in
# the functions found in cylind_simulation_functions.R 
# At the end of the scripts there is a suite of tests for verifying the accuracy 
# of the convolutions and the approximation; by default these are not run.
# IN: N/A
# OUT: ../ValidationData/conv_test_simulations.RData
# TEST: Does not overwrite existing data

test_mode=T
library(tidyverse)
library(BayesianOHC)

rangefun_stat=function(x) rep(x,nlocs)
rangefun_nonstat=function(x) (cos(locs)+1)*x/2+x/2

locs=sort(runif(100,-pi,pi))
nranges=100
nrep=100
to_effrange<-function(x) sqrt(-log(0.05)*x)

#Sequence up to 120 degree effective range
true_range_seq=seq(.001,deg2rad(120),length=nranges)^2/(-log(0.05))

stationary_simulations=
  cylind_approx_simulation(true_range_seq,rangefun_stat,nrep,locs)
nonstationary_simulations=
  cylind_approx_simulation(true_range_seq,rangefun_nonstat,nrep,locs)

if(!test_mode){
  save(stationary_simulations,nonstationary_simulations,
       file='../ValidationData/conv_test_simulations.RData')
}

run_tests=F
#Below are tests to ensure that the convolutions and approximation are
#giving the expected results
if(run_tests){
  ###Test that euclidean_correlation_convolution is equal to numeric convolution
  for(i in 1:10){
    kernel_euc<-function(u,x,range){
      return(exp(-(u-x)^2/range))
    }
    convfun_euc<-function(x,y,range1,range2,low,up){
      myseq=seq(low,up,length=1000000)
      return(mean(kernel_euc(myseq,x,range1)*
                    kernel_euc(myseq,y,range2))*(up-low))
    }
    x=runif(1,-pi,pi)
    y=x+runif(1,-pi,pi)
    range1=runif(1,0,4)
    range2=runif(1,0,7)
    lower=-10
    upper=10
    
    val1=convfun_euc(x,y,range1,range2,lower,upper)
    scaling_factor=sqrt(pi/2)*sqrt(sqrt(range1*range2))
    val2=scaling_factor*euclidean_correlation_convolution(x,y,2*range1,2*range2)
    
    print(val1/val2-1)
  }
  
  ###Test that cylindrical_correlation_exact is equal to actual convolution
  for(i in 1:10){
    kernel_gc<-function(u,x,range2){
      return(exp(-mygc_dist(u,x)^2/(range2)))
    }
    convfun_gc<-function(x,y,range1,range2,low,up){
      myseq=seq(low,up,length=1000000)
      return(mean(kernel_gc(myseq,x,range1)*kernel_gc(myseq,y,range2))*(up-low))
    }
    
    x=runif(1,-pi,pi)
    y=runif(1,-pi,pi)
    range1=runif(1,0,4)
    range2=runif(1,0,7)
    
    val1=convfun_gc(x,y,range1,range2,-pi,pi)
    scaling_term=(sqrt(pi/2)*sqrt(sqrt(range1*range2)))
    val2=scaling_term*cylindrical_correlation_exact(x,y,2*range1,2*range2)
    
    print(val1/val2-1)
  }
  
  ###Check that cylindrical_correlation_exact approximates cylindrical_correlation_gaussian
  ###Note, not all will  be near zero, but most should be
  for(i in 1:10){
    x=runif(1,-pi,pi)
    y=runif(1,-pi,pi)
    range1=runif(1,0,1)
    range2=runif(1,0,2)
    val1=cylindrical_correlation_gaussian(x,y,range1,range2)
    val2=cylindrical_correlation_exact(x,y,range1,range2)
    print(val1/val2-1)
  }
  
}


