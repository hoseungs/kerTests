library(mvtnorm)

simul_power1 = function(type, m, n, d, a, b, permut, iter, alpha) {
  # new test
  if (type==1) {
    res_fGPK_appr = res_fGPK_perm = res_fGPKM_appr = res_fGPKM_perm = res_GPK_perm = 0
    for (i in 1:iter) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = b*A)
      sigma = med_sigma(X,Y)
      res = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
      if (res$pval$fGPK_appr <= alpha) {
        res_fGPK_appr = res_fGPK_appr + 1
      }
      if (res$pval$fGPK_perm <= alpha) {
        res_fGPK_perm = res_fGPK_perm + 1
      }
      if (res$pval$fGPKM_appr <= alpha) {
        res_fGPKM_appr = res_fGPKM_appr + 1
      }
      if (res$pval$fGPKM_perm <= alpha) {
        res_fGPKM_perm = res_fGPKM_perm + 1
      }
      if (res$pval$GPK_perm <= alpha) {
        res_GPK_perm = res_GPK_perm + 1
      }
    }
    res_power = list()
    res_power$fGPK_appr = res_fGPK_appr
    res_power$fGPK_perm = res_fGPK_perm
    res_power$fGPKM_appr = res_fGPKM_appr
    res_power$fGPKM_perm = res_fGPKM_perm
    res_power$GPK_perm = res_GPK_perm
  }
  # MMD Bootstrap
  if (type==2) {
    res_power = 0
    for (i in 1:iter) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = b*A)
      sigma = med_sigma(X,Y)
      res = mmd(X, Y, sigma, B=permut)
      if (res <= alpha) {
        res_power = res_power + 1
      }
    }
  }
  # AG
  if (type==3) {
    res_power = 0
    for (i in 1:iter) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = b*A)
      N = n+m
      Nresample = 1000
      Permutation = t(replicate(Nresample,sample(1:N,N,replace=F)))
      d2mat = (as.matrix(dist(rbind(X,Y),method="euclidean",diag=TRUE,upper=TRUE)))^2
      upper_bound = min(2/d*log(n),4);
      lower_bound = -2
      svector = exp(seq(lower_bound,upper_bound,by=1/4))
      U = length(svector)
      
      Var_vector = NULL
      for (s in svector){
        Var_vector = c(Var_vector,Var_ts(d2mat,s))
        Var_vector = c(Var_vector,max(Var_ts(d2mat,s),1/n^2))
      }
      nnld_stat_vector = sapply(svector,FUN=nnld_MMD2_ts,d2mat=d2mat,permutation=1:N,n=n)
      snld_stat_vector = nnld_stat_vector/sqrt(2*Var_vector)/(1/n+1/m)
      snld = max(snld_stat_vector)
      
      snld_re = rep(0,Nresample)
      
      for (r in 1:Nresample)
      {
        permutation = Permutation[r,]
        nnld_stat_vector_re = sapply(svector,FUN=nnld_MMD2_ts,d2mat=d2mat,permutation=permutation,n=n)
        snld_stat_vector_re = nnld_stat_vector_re/sqrt(2*Var_vector)/(1/n+1/m)
        snld_re[r] = max(snld_stat_vector_re)
      }
      res = (sum(snld_re>snld)+1)/(Nresample+1)
      if (res <= alpha) {
        res_power = res_power + 1
      }
    }
  }
  # ND1
  if (type==4) {
    res_power = 0
    for (i in 1:iter) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = b*A)
      v = (m-1)*(n-1) + m*(n-3)/2 + n*(n-3)/2
      res = Test_stat(X,Y,"G-induced",1)
      res = pt(-res, df=v)
      if (res <= alpha) {
        res_power = res_power + 1
      }
    }
  }
  # ND2
  if (type==5) {
    res_power = 0
    for (i in 1:iter) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = b*A)
      v = (m-1)*(n-1) + m*(n-3)/2 + n*(n-3)/2
      res = Test_stat(X,Y,"G-induced",2)
      res = pt(-res, df=v)
      if (res <= alpha) {
        res_power = res_power + 1
      }
    }
  }
  # TM
  if (type==6) {
    res_power = 0
    for (i in 1:iter) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = b*A)
      res = compute.T.k(X, Y, kernel="G")
      res = pnorm(-res$T.k)
      if (res <= alpha) {
        res_power = res_power + 1
      }
    }
  }
  
  return(res_power)
}
