kertests = function(X, Y, sigma, r1=1.2, r2=0.8, perm=0) {
  m = nrow(X)
  n = nrow(Y)

  Z = rbind(X,Y)
  N = nrow(Z)

  K = exp(-as.matrix(dist(Z)^2)/sigma^2/2) - diag(N) # kernel matrix

  Kx = sum(K[1:m,1:m])/m/(m-1)
  Ky = sum(K[(m+1):N,(m+1):N])/n/(n-1)
  Kxy = sum(K[1:m,(m+1):N])/m/n

  mu_Kx = sum(K)/N/(N-1)
  mu_Ky = sum(K)/N/(N-1)
  mu_Kxy = sum(K)/N/(N-1)

  A = sum(K^2)
  B = sum(rowSums(K)^2) - A
  C = sum(K)^2 - 2*A - 4*B

  p1 = m*(m-1)/N/(N-1)
  p2 = p1*(m-2)/(N-2)
  p3 = p2*(m-3)/(N-3)

  q1 = n*(n-1)/N/(N-1)
  q2 = q1*(n-2)/(N-2)
  q3 = q2*(n-3)/(N-3)

  var_Kx = (2*A*p1 + 4*B*p2 + C*p3)/m/m/(m-1)/(m-1) - mu_Kx^2
  var_Ky = (2*A*q1 + 4*B*q2 + C*q3)/n/n/(n-1)/(n-1) - mu_Ky^2
  cov_Kx_Ky = C/N/(N-1)/(N-2)/(N-3) - mu_Kx*mu_Ky

  # test statistic GPK
  COV = matrix(c(var_Kx,cov_Kx_Ky,cov_Kx_Ky,var_Ky), nrow=2)
  Sinv = solve(COV)
  kmv = c(Kx-mu_Kx, Ky-mu_Ky)
  GPK = as.numeric(kmv%*%Sinv%*%kmv)

  # test statistic Z_D
  u.D = m*(m-1)
  v.D = -n*(n-1)
  mean.D = mu_Kx*u.D + mu_Ky*v.D
  var.D = (u.D^2)*var_Kx + (v.D^2)*var_Ky + 2*u.D*v.D*cov_Kx_Ky
  Z.D = (Kx*u.D + Ky*v.D - mean.D)/sqrt(var.D)

  # test statistic Z_W1
  u.W1 = r1*m/N
  v.W1 = n/N
  mean.W1 = mu_Kx*u.W1 + mu_Ky*v.W1
  var.W1 = var_Kx*u.W1^2 + var_Ky*v.W1^2 + 2*u.W1*v.W1*cov_Kx_Ky
  Z.W1 = (Kx*u.W1 + Ky*v.W1 - mean.W1)/sqrt(var.W1)

  # test statistic Z_W2
  u.W2 = r2*m/N
  v.W2 = n/N
  mean.W2 = mu_Kx*u.W2 + mu_Ky*v.W2
  var.W2 = var_Kx*u.W2^2 + var_Ky*v.W2^2 + 2*u.W2*v.W2*cov_Kx_Ky
  Z.W2 = (Kx*u.W2 + Ky*v.W2 - mean.W2)/sqrt(var.W2)
  
  temp_approx = sort( c(pnorm(-Z.W1), pnorm(-Z.W2), 2*pnorm(-abs(Z.D))) )
  fGPK_appr = 3*min( temp_approx[1], temp_approx[2], temp_approx[3] )
  
  temp_approx1 = sort( c(pnorm(-Z.W1), pnorm(-Z.W2)) )
  fGPKM_appr = 2*min( temp_approx1[1], temp_approx1[2] )

  fGPK_Simes_appr = min( 3*temp_approx[1], 1.5*temp_approx[2], temp_approx[3] )
  fGPKM_Simes_appr = min( 2*temp_approx1[1], temp_approx1[2] )

  result = list()
  result$teststat$GPK = GPK
  result$teststat$ZW1 = Z.W1
  result$teststat$ZW2 = Z.W2
  result$teststat$ZD = Z.D
  result$pval$fGPK_appr = min(1,fGPK_appr)
  result$pval$fGPKM_appr = min(1,fGPKM_appr)
  result$pval$fGPK_Simes_appr = min(1,fGPK_Simes_appr)
  result$pval$fGPKM_Simes_appr = min(1,fGPKM_Simes_appr)


  if (perm>0) {
    temp1 = temp2 = temp3 = temp4 = rep(0, perm)
    for (i in 1:perm) {
      id = sample(N, replace = FALSE)
      K_i = K[id,id]

      Kx_i = sum(K_i[1:m,1:m])/m/(m-1)
      Ky_i = sum(K_i[(m+1):N,(m+1):N])/n/(n-1)
      Kxy_i = sum(K_i[1:m,(m+1):N])/m/n

      kmv_i = c(Kx_i-mu_Kx, Ky_i-mu_Ky)
      GPK_i = as.numeric(kmv_i%*%Sinv%*%kmv_i)

      Z.D_i = (Kx_i*u.D + Ky_i*v.D - mean.D)/sqrt(var.D)
      Z.W1_i = (Kx_i*u.W1 + Ky_i*v.W1 - mean.W1)/sqrt(var.W1)
      Z.W2_i = (Kx_i*u.W2 + Ky_i*v.W2 - mean.W2)/sqrt(var.W2)

      temp1[i] = GPK_i
      temp2[i] = Z.D_i
      temp3[i] = Z.W1_i
      temp4[i] = Z.W2_i
    }
    GPK_perm = length(which(temp1>=GPK))/perm

    perm_pval_Z.D = 2*length(which(temp2>=abs(Z.D)))/perm
    perm_pval_Z.W1 = length(which(temp3>=Z.W1))/perm
    perm_pval_Z.W2 = length(which(temp4>=Z.W2))/perm
    
    temp_perm = sort(c(perm_pval_Z.W1, perm_pval_Z.W2, perm_pval_Z.D))
    fGPK_perm = 3*min( temp_perm[1], temp_perm[2], temp_perm[3] )
    
    temp_perm1 = sort(c(perm_pval_Z.W1, perm_pval_Z.W2))
    fGPKM_perm = 2*min( temp_perm1[1], temp_perm1[2])

    fGPK_Simes_perm = min( 3*temp_perm[1], 1.5*temp_perm[2], temp_perm[3] )
    fGPKM_Simes_perm = min( 2*temp_perm1[1], temp_perm1[2])

    result$pval$GPK_perm = min(1,GPK_perm)
    result$pval$fGPK_perm = min(1,fGPK_perm)
    result$pval$fGPKM_perm = min(1,fGPKM_perm)
    result$pval$fGPK_Simes_perm = min(1,fGPK_Simes_perm)
    result$pval$fGPKM_Simes_perm = min(1,fGPKM_Simes_perm)
  }

  return(result)
}

