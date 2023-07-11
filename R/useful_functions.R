# function to get the median of all pairwise distances among observations for the bandwidth
med_sigma = function(X, Y) {
  aggre = rbind(X,Y)
  med = median(dist(aggre)^2)
  return(sqrt(med*0.5))
}

# covariance structure in data
Sig = function(val, d) {
  cov = matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      cov[i,j] = val^(abs(i-j))
    }
  }
  return(cov)
}


# MMD-Bootstrap
mmd = function(X, Y, sigma, B=0) {
  m = nrow(X)
  n = nrow(Y)
  
  Z = rbind(X,Y)
  N = nrow(Z)
  
  K = exp(-as.matrix(dist(Z)^2)/sigma^2/2) - diag(N)
  
  Kx = sum(K[1:m,1:m])/m/(m-1)
  Ky = sum(K[(m+1):N,(m+1):N])/n/(n-1)
  Kxy = sum(K[1:m,(m+1):N])/m/n
  
  MMD = Kx + Ky - 2*Kxy
  
  if (B>0) {
    temp =  rep(0, B)
    for (i in 1:B) {
      id = sample(N, replace = TRUE)
      K_i = K[id,id]
      
      Kx_i = sum(K_i[1:m,1:m])/m/(m-1)
      Ky_i = sum(K_i[(m+1):N,(m+1):N])/n/(n-1)
      Kxy_i = sum(K_i[1:m,(m+1):N])/m/n
      
      MMD_i = Kx_i + Ky_i - 2*Kxy_i
      
      temp[i] = MMD_i
    }
    pval = length(which(temp>=MMD))/B
  }
  
  return(pval)
}


# MMD-Pearson
library(PearsonDS)

mPearson = function(X, Y, sigma) {
  m = nrow(X)
  
  Z = rbind(X,Y)
  N = 2*m
  
  K = exp(-as.matrix(dist(Z)^2)/sigma^2/2)
  
  Kx = K[1:m,1:m]
  Ky = K[(m+1):N,(m+1):N]
  Kxy = K[1:m,(m+1):N]
  
  MMDf = Kx + Ky - Kxy - t(Kxy)
  MMDf = MMDf - diag(diag(MMDf))
  
  test.stat = sum(MMDf)/(m-1)
  
  # 2nd moment
  MMDg = MMDf^2
  MMD = MMDg/m/(m-1)
  preVarMMD = sum(MMD)
  U_m2 = preVarMMD*2/m/(m-1)
  
  # 3rd moment
  pre3rdMMD = 0
  for (i1 in 1:m) {
    for (i2 in 1:m) {
      if (i1 != i2) {
        intMean = 0
        for (i3 in 1:m) {
          if (i3 != i1 && i3 != i2) {
            intMean = intMean + MMDf[i1,i3]*MMDf[i2,i3]
          }
        }
        intmean = intMean/(m-2)
        pre3rdMMD = pre3rdMMD + intmean*MMDf[i1,i2]
      }
    }
  }
  pre3rdMMD = pre3rdMMD/m/(m-1)
  U_m3 = 8*(m-2)/(m*(m-1))^2 * pre3rdMMD
  
  # Rescale 2nd and 3rd moments, since they were computed for MMD and
  # must be applied to m*MMD
  U_m3 = m^3 * U_m3
  U_m2 = m^2 * U_m2
  
  # Approx 4th moment
  # kurtosis >= skewness^2 + 1
  kurt = 2*((U_m3)^2/U_m2^3 + 1)
  
  # Sample the distribution
  U_moments = c(0, sqrt(U_m2), U_m3/U_m2^(3/2), kurt)
  
  samp = rpearson(100, moments = U_moments)
  pval = length(which(samp>test.stat))/100
  
  return( pval )
}

