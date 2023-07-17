library(mvtnorm)


### $ d_i = d_0 $ for i=1,..,p (the group sizes; see Section 3 in Chakraborty & Zhang (2021)) #########
### Here we allow d_0 = 1 or 2. In Examples 6.1 & 6.2 in the paper, we use d_0 = 1. In Example 6.3, we use d_0 = 2.


## gam2 : computes the median distance (excluding zeroes) among rows of the aggregated data matrix 

## Input : 
## u : a numeric N*p data matrix (N=n+m= aggregated sample size, p= dimension of the random vector).

gam2 <- function(u){  
  D = as.matrix(dist(u,"euclidean")); D_vec = as.vector(D)
  val = median(D_vec[D_vec!=0])
  return(val)
}


## gam_grp : computes the groupwise median distance (excluding zeroes) among rows of the aggregated
## data matrix

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector) for X.
## y : a numeric m*p data matrix (m=sample size, p= dimension of the random vector) for Y.
## d_0 : dimension of the groups (1 or 2)

gam_grp <- function(x,y,d_0){  
  comb=as.matrix(rbind(x,y))
  combt=t(comb)
  if(d_0==1){
    gnew = as.vector(apply(combt, 1, gam2))
  }
  if(d_0==2){
    gnew = as.vector(unlist(lapply(lapply(split.data.frame(combt,rep(1:(dim(combt)[1]/2),each=2)), t), gam2))) 
  }
  return(gnew)
}


## define some functions

f_fun <- function(x){
  return(sqrt(sum(x^2)))
}

f_list <- function(x){
  lp <- length(x)
  li=split(x,rep(1:(lp/2),each=2))
  v_vec <- as.vector(unlist(lapply(li,f_fun)))
  return(sqrt(sum(v_vec)))
}


## XXdist : computes the distance matrix for X, given one of the 6 distance metrics
## specified in Section 6.1, Chakraborty & Zhang (2021).

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector) for X.
## z : a numeric m*p data matrix (m=sample size, p= dimension of the random vector) for Y (required to compute
## the tuning/bandwidth parameter for the aggregated sample according to the median heuristic).
## type : either of "K", "L-induced", "G-induced", "E", "MMD-E" and "MMD-G", to choose the 
## metric from among the 6 listed metrics in Section 6.1, Chakraborty & Zhang (2021).
## d_0 : dimension of the groups (1 or 2)

XXdist <- function(x,z,type,d_0){   
  x1 = t(x)                    
  n = dim(x)[1]; p=dim(x)[2]
  d = matrix(NA, n, n)
  
  if(type == "E") { for(i in 1:n) d[i,] = sqrt( apply((x1 - x1[,i])^2, 2, sum) ) }  ## for usual Euclidean ED
  
  if(type == "K") { 
    
    if(d_0==1){
      for(i in 1:n) d[i,] = sqrt( apply(abs(x1 - x1[,i]), 2, sum) )
    }
    if(d_0==2){
      for(i in 1:n) d[i,] = apply(x1 - x1[,i], 2, f_list) 
    }
  }
  ## for \cal{E} with d_i=1 and \rho_i is the Euclidean distance
  
  g=gam2(rbind(x,z))
  
  if(type == "MMD-E")     ## for MMD with Laplace kernel
  { 
    for(i in 1:n) d[i,] = sqrt( apply((x1 - x1[,i])^2, 2, sum) ) 
    d = exp(-d/g) 
    diag(d) = 0
  }
  if(type == "MMD-G")     ## for MMD with Gaussian kernel
  { 
    for(i in 1:n) d[i,] = sqrt( apply((x1 - x1[,i])^2, 2, sum) ) 
    d = exp(-d^2/g^2/2)
    diag(d) = 0
  }
  
  gnew=gam_grp(x,z,d_0)
  
  if(type=="L-induced"){  ## for \cal{E} with d_i=1 and \rho_i is the distance induced by the Laplace kernel
    for(i in 1:n){
      #if(i%%2==0) print(i)
      d[i,i]=0
      if(i<n){
        for(l in (i+1):n){
          
          if(d_0==1){
            d[i,l]=sqrt( sum( 2-2*exp(-abs(x1[,i] - x1[,l])/gnew )))
          }
          if(d_0==2){
            d[i,l]=sqrt( sum( 2-2*exp(- (as.vector(unlist(lapply(split(x1[,i] - x1[,l],rep(1:(p/2),each=2)),f_fun))))/gnew ))) 
          }
          
          d[l,i]=d[i,l]
        }
      }
    }
  }
  
  if(type=="G-induced"){  ## for \cal{E} with d_i=1 and \rho_i is the distance induced by the Gaussian kernel
    for(i in 1:n){
      #if(i%%10==0) print(i)
      d[i,i]=0
      if(i<n){
        for(l in (i+1):n){
          
          if(d_0==1){
            d[i,l]=sqrt( sum( 2-2*exp(-(x1[,i] - x1[,l])^2 /2/(gnew)^2 )))
          }
          if(d_0==2){
            d[i,l]=sqrt( sum( 2-2*exp(-(as.vector(unlist(lapply(split(x1[,i] - x1[,l],rep(1:(p/2),each=2)),f_fun))))^2 /2/(gnew)^2 )))
          }
          
          d[l,i]=d[i,l]
        }
      }
    }
  }
  
  return(d)
}


## XYdist : computes the distance matrix between X & Y, given one of the 6 distance metrics
## specified in Section 6.1, Chakraborty & Zhang (2021).

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector) for X.
## y : a numeric m*p data matrix (m=sample size, p= dimension of the random vector) for Y
## type : either of "K", "L-induced", "G-induced", "E", "MMD-E" and "MMD-G", to choose the 
## metric from among the 6 listed metrics in Section 6.1, Chakraborty & Zhang (2021).
## d_0 : dimension of the groups (1 or 2).

XYdist <- function(x,y,type,d_0){    
  
  x1 <- t(x) ; y1 <- t(y); p=dim(x)[2]
  n1 <- dim(x)[1] ; n2 <- dim(y)[1]
  dxy <- matrix(NA, n1, n2); dxy1 <- matrix(NA, n1, n2)
  
  if(type == "E") { for(i in 1:n1) dxy[i,] = sqrt( apply((y1 - x1[,i])^2, 2, sum) ) }
  if(type == "K") {
    
    if(d_0==1){
      for(i in 1:n1) dxy[i,] = sqrt( apply(abs(y1 - x1[,i]), 2, sum) )
    }
    if(d_0==2){
      for(i in 1:n1) dxy[i,] = apply(y1 - x1[,i], 2, f_list)
    }
  }
  
  g=gam2(rbind(x,y))
  
  if(type == "MMD-E")     ## for MMD with Laplace kernel
  { 
    for(i in 1:n1) dxy1[i,] = sqrt( apply((y1 - x1[,i])^2, 2, sum) ) 
    dxy = exp(-dxy1/g) 
  }
  if(type == "MMD-G")     ## for MMD with Gaussian kernel
  { 
    for(i in 1:n1) dxy1[i,] = sqrt( apply((y1 - x1[,i])^2, 2, sum) ) 
    dxy = exp(-dxy1^2/g^2/2)
  }
  
  gnew1=gam_grp(x,y,d_0)
  
  if(type=="L-induced"){ ## for \cal{E} with d_i=1 and \rho_i is the distance induced by the Laplace kernel
    for(i in 1:n1){
      for(l in 1:n2){
        
        if(d_0==1){
          dxy[i,l]=sqrt( sum( 2-2*exp(-abs(x1[,i] - y1[,l])/gnew1 )))
        }
        if(d_0==2){
          dxy[i,l]=sqrt( sum( 2-2*exp(-(abs(as.vector(unlist(lapply(split(x1[,i] - y1[,l],rep(1:(p/2),each=2)),f_fun)))))/gnew1 )))
        }
        
      }
    }
  }
  
  if(type=="G-induced"){ ## for \cal{E} with d_i=1 and \rho_i is the distance induced by the Gaussian kernel
    for(i in 1:n1){
      for(l in 1:n2){
        
        if(d_0==1){
          dxy[i,l]=sqrt( sum( 2-2*exp(-(x1[,i] - y1[,l])^2 /2/(gnew1)^2 )))
        }
        if(d_0==2){
          dxy[i,l]=sqrt( sum( 2-2*exp(-(as.vector(unlist(lapply(split(x1[,i] - y1[,l],rep(1:(p/2),each=2)),f_fun))))^2 /2/(gnew1)^2 )))
        }
        
      }
    }
  }
  
  return(dxy)
}


## ED : computes the U-statistic type estimator of E(X,Y) defined in Section 4 in 
## Chakraborty and Zhang (2021).

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector) for X.
## y : a numeric m*p data matrix (m=sample size, p= dimension of the random vector) for Y
## type : either of "K", "L-induced", "G-induced", "E", "MMD-E" and "MMD-G", to choose the 
## metric from among the 6 listed metrics in Section 6.1, Chakraborty & Zhang (2021).
## d_0 : dimension of the groups (1 or 2).

ED <- function(x,y,type,d_0){    
  n1 <- dim(x)[1] ; n2 <- dim(y)[1]
  
  Axx=XXdist(x,y,type,d_0)
  Ayy=XXdist(y,x,type,d_0)
  Axy=XYdist(x,y,type,d_0)
  
  if(type=="K" || type=="E" || type=="L-induced" || type=="G-induced"){
    ed <- sum(Axy) * (2/(n1*n2)) - sum(Axx) *(1/(n1*(n1-1))) - sum(Ayy) *(1/(n2*(n2-1)))
  }
  else{
    ed <- sum(Axx) *(1/(n1*(n1-1))) + sum(Ayy) *(1/(n2*(n2-1))) - sum(Axy) * (2/(n1*n2)) 
  }
  l=list(ed,Axx,Ayy,Axy)
  return(l)
}

## cdCov : computes the quantity `Cross Distance Covariance' defined in Section 4 in Chakraborty and Zhang (2021)
## given the distance matrix between X & Y

## Input :
## D : a numeric distance matrix between X & Y

cdCov <- function(D){  
  n1=dim(D)[1]; n2=dim(D)[2]
  #D <- XYdist(x,y,type)
  R <- rowMeans(D) ; C <- colMeans(D) ; T <- mean(D)
  Rm <- matrix(rep(R,n2),n1,n2)
  Cm <- matrix(rep(C,n1),n1,n2, byrow=TRUE)
  Tm <- matrix(rep(T, n1*n2), n1, n2)
  Dhat <- D - Rm - Cm + Tm
  A <- sum(Dhat*Dhat)/((n1-1)*(n2-1))
  return(A)
}

## u.center : computes the U-centered version of a given distance matrix, given in equation (1.1) in 
## Chakraborty and Zhang (2021).

## Input : 
## A : a numeric distance matrix.

u.center = function(A){    
  n = dim(A)[1]                 
  R = rowSums(A)
  C = colSums(A)
  T = sum(A)
  r = matrix(rep(R,n),n,n)/(n-2)
  c = t(matrix(rep(C,n),n,n))/(n-2)
  t = matrix(T/(n-1)/(n-2),n,n)
  UA = -(A-r-c+t)
  diag(UA)=0
  return(UA)
}

## MdCov.U : computes the U-statistic type estimator of D(X,X) defined in Section 5 in 
## Chakraborty and Zhang (2021).

## Input : 
## B : a numeric distance matrix for the samples of X.

MdCov.U = function(B) {   
  n = dim(B)[1] 
  A1 = u.center(B)
  return( sum(A1*A1)/n/(n-3) )
}



###### -------------------------- Main function ------------------------------- ######


## Test_stat : computes the value of the test statistic $T_{n,m}$ defined 
## in Section 4.1 in  Chakraborty and Zhang (2021).

## Input : 
## x : a numeric n*p data matrix (n=sample size, p= dimension of the random vector) for X.
## y :  a numeric m*p data matrix (m=sample size, p= dimension of the random vector) for Y.
## type : either of "K", "L-induced", "G-induced", "E", "MMD-E" and "MMD-G", to choose the 
## metric from among the 6 listed metrics in Section 6.1, Chakraborty & Zhang (2021).
## d_0 : dimension of the groups (1 or 2).

Test_stat <- function(x,y,type,d_0){   ## computes the test statistic, inputs sample matrices & metric type
  n1 <- dim(x)[1] ; n2 <- dim(y)[1]
  L=ED(x,y,type,d_0)
  Num <- L[[1]]
  AXX=L[[2]]; AYY=L[[3]]; AXY=L[[4]]
  scalar.den <- 1/(n1*n2) + 1/(2*n1*(n1-1)) + 1/(2*n2*(n2-1))
  
  numSnm <- 4*(n1-1)*(n2-1)*cdCov(AXY) + 4*n1*(n1-3)/2 * MdCov.U(AXX) + 4*n2*(n2-3)/2 * MdCov.U(AYY)
  denSnm <- (n1-1)*(n2-1) + n1*(n1-3)/2 + n2*(n2-3)/2
  Snm = numSnm/denSnm
  
  Den <- sqrt(scalar.den * Snm)
  return(Num/Den)
}
