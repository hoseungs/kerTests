library(mvtnorm)

simul1 = function(type, m, n, d, mu, adj, permut) {
  # new test
  if (type==1) {
    A = Sig(0.4,d)
    X = rmvnorm(m, mean = rep(0, d), sigma = A)
    Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
    sigma = med_sigma(X,Y)
    a = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
  }
  # MMD Bootstrap
  if (type==2) {
    A = Sig(0.4,d)
    X = rmvnorm(m, mean = rep(0, d), sigma = A)
    Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
    sigma = med_sigma(X,Y)
    a = mmd(X, Y, sigma, B=permut)
  }
  # AG
  if (type==3) {
    A = Sig(0.4,d)
    X = rmvnorm(m, mean = rep(0, d), sigma = A)
    Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
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
    a = (sum(snld_re>snld)+1)/(Nresample+1)
  }
  # ND1
  if (type==4) {
    A = Sig(0.4,d)
    X = rmvnorm(m, mean = rep(0, d), sigma = A)
    Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
    v = (m-1)*(n-1) + m*(n-3)/2 + n*(n-3)/2
    a = Test_stat(X,Y,"G-induced",1)
    a = pt(-a, df=v)
  }
  # ND2
  if (type==5) {
    A = Sig(0.4,d)
    X = rmvnorm(m, mean = rep(0, d), sigma = A)
    Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
    v = (m-1)*(n-1) + m*(n-3)/2 + n*(n-3)/2
    a = Test_stat(X,Y,"G-induced",2)
    a = pt(-a, df=v)
  }
  # TM
  if (type==6) {
    A = Sig(0.4,d)
    X = rmvnorm(m, mean = rep(0, d), sigma = A)
    Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
    a = compute.T.k(X, Y, kernel="G")
    a = pnorm(-a$T.k)
  }
  
  return(a)
}
