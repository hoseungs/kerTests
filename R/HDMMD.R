
##################################################################
#####                                                        #####
#####       TWO-SAMPLE HIGH-DIMENSIONAL MMD-BASED TEST       #####  
#####                                                        #####
#####       AUTHORS: HANJIA GAO, XIAOFENG SHAO               #####
#####       INSTITUTE: UNIVERSITY OF ILLINOIS                #####
#####                  AT URBANA-CHAMPAIGN                   #####                    
##################################################################

## Given independent random samples X=(X_1,...,X_n) and Y=(Y_1,...,Y_m)
## Let Z=(X,Y) denote the pooled sample over X and Y
## Test the null hypothesis H_0: X =^d Y
## The test statistic T_{n,m,p}^{k}(X,Y) is defined as Equ(3.6) in the manuscript
## Reject the null at level alpha if T_{n,m,p}^{k}(X,Y) > qnorm(1-alpha)


compute.gamma = function(Z,kernel){
  switch (kernel,
          G = median(stats::dist(Z, method = "euclidean")^2),
          L = median(stats::dist(Z, method = "euclidean"))
  )
}


compute.a0 = function(kernel){
  switch (kernel,
          L2 = 0,
          G = -1,
          L = -1
  )
}


compute.a.st.k = function(Z,kernel,k.gamma=1){
  switch(kernel,
         L2 = as.matrix(stats::dist(Z, method = "euclidean")),
         L1 = as.matrix(stats::dist(Z, method = "manhattan")),
         G = -exp(-as.matrix(stats::dist(Z, method = "euclidean"))^2/(2*k.gamma)),
         L = -exp(-as.matrix(stats::dist(Z, method = "euclidean"))/k.gamma)
  )
}


u.centered.dist = function(a.st.k){
  n = nrow(a.st.k)
  A.1 = a.st.k
  A.2 = kronecker(cbind(rep(1,n)), rbind(colSums(a.st.k)))/(n-2)
  A.3 = kronecker(cbind(rowSums(a.st.k)), rbind(rep(1,n)))/(n-2)
  A.4 = matrix(sum(a.st.k), nrow = nrow(a.st.k), ncol = ncol(a.st.k))/((n-1)*(n-2))
  return(A.1 - A.2 - A.3 + A.4)
}


compute.v.k = function(a.st.k){
  A.st.k = u.centered.dist(a.st.k)
  n = nrow(A.st.k)
  return((sum(A.st.k^2)-sum(diag(A.st.k^2)))/(n*(n-3)))
}


compute.e.k = function(a.st.k,n,m){
  c = matrix(1/(n*m), n+m, n+m)
  c[1:n,1:n] = -1/(n*(n-1))
  c[(n+1):(n+m),(n+1):(n+m)] = -1/(m*(m-1))
  diag(c) = 0
  return(sum(c*a.st.k))
}


## returns the test statistic T_{n,m,p}^{k}(X,Y)
## the output is a list containing 
##   (1) the sample estimate V_{n,m}^{k\ast}(X,Y) **WITH bias-correction**, named by v.k
##   (2) the sample estimate E_{n,m}^{k}(X,Y), named by e.k
##   (3) the test statistic T_{n,m,p}^{k}(X,Y), named by T.k
##   (4) the computational time counted in secs, named by time.k
## rejects the null at level alpha if T_{n,m,p}^{k}(X,Y) > qnorm(1-alpha)
compute.T.k = function(X,Y,kernel){
  st = Sys.time()
  n = nrow(X)
  m = nrow(Y)
  rho = n / (n+m)
  Z = rbind(X,Y)
  
  if ((kernel == "G") || (kernel == "L")){
    k.gamma = compute.gamma(Z,kernel)
  }
  
  a.st.k = compute.a.st.k(Z,kernel,k.gamma)
  v.k.Z = compute.v.k(a.st.k) - compute.a0(kernel)^2/(n+m-1)/(n+m-3)
  c.nm = 2/(n*(n-1)) + 4/(n*m) + 2/(m*(m-1))
  e.k = compute.e.k(a.st.k,n,m)
  ed = Sys.time()
  
  return(list(v.k = v.k.Z,
              e.k = e.k,
              T.k = e.k/sqrt(c.nm * v.k.Z),
              time.k = as.numeric(difftime(ed,st),units="secs")))
}

