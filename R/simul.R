library(mvtnorm)
library(ade4)
library(Ball)
library(gTests)
library(Ecume)
library(expm)

simul = function(type, m, n, d, dat, mu, adj, permut) {
  # new test
  if (type==1) {
    # Gaussian
    if (dat==1) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A) 
      Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
      sigma = med_sigma(X,Y)
      a = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
    }
    # t
    if (dat==2) {
      A = Sig(0.4,d)
      B = Sig(0.4,d)
      X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
      Y = rmvt(n, sigma = adj*A, delta = rep(mu, d), df = 20)
      sigma = med_sigma(X,Y)
      a = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
    }
    # chi
    if (dat==3) {
      A = sqrtm( Sig(0.4,d) )
      B = sqrtm( adj*Sig(0.4,d) )
      X = t( A%*%matrix(rchisq(d*m, df=3), d) )
      Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(mu, d, n) )
      sigma = med_sigma(X,Y)
      a = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
    }
  }
  # MMD Bootstrap
  if (type==2) {
    if (dat==1) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
      sigma = med_sigma(X,Y)
      a = mmd(X, Y, sigma, B=permut)
    }
    if (dat==2) {
      A = Sig(0.4,d)
      B = Sig(0.4,d)
      X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
      Y = rmvt(n, sigma = adj*A, delta = rep(mu, d), df = 20)
      sigma = med_sigma(X,Y)
      a = mmd(X, Y, sigma, B=permut)
    }
    if (dat==3) {
      A = sqrtm( Sig(0.4,d) )
      B = sqrtm( adj*Sig(0.4,d) )
      X = t( A%*%matrix(rchisq(d*m, df=3), d) )
      Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(mu, d, n) )
      sigma = med_sigma(X,Y)
      a = mmd(X, Y, sigma, B=permut)
    }
  }
  # MMD Pearson
  if (type==3) {
    if (dat==1) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
      sigma = med_sigma(X,Y)
      a = mPearson(X, Y, sigma)
    }
    if (dat==2) {
      A = Sig(0.4,d)
      B = Sig(0.4,d)
      X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
      Y = rmvt(n, sigma = adj*A, delta = rep(mu, d), df = 20)
      sigma = med_sigma(X,Y)
      a = mPearson(X, Y, sigma)
    }
    if (dat==3) {
      A = sqrtm( Sig(0.4,d) )
      B = sqrtm( adj*Sig(0.4,d) )
      X = t( A%*%matrix(rchisq(d*m, df=3), d) )
      Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(mu, d, n) )
      sigma = med_sigma(X,Y)
      a = mPearson(X, Y, sigma)
    }
  }
  # BT
  if (type==4) {
      if (dat==1) {
          A = Sig(0.4,d)
          X = rmvnorm(m, mean = rep(0, d), sigma = A)
          Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
          a = bd.test(X, Y, method = "limit")
      }
      if (dat==2) {
          A = Sig(0.4,d)
          B = Sig(0.4,d)
          X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
          Y = rmvt(n, sigma = adj*A, delta = rep(mu, d), df = 20)
          a = bd.test(X, Y, method = "limit")
      }
      if (dat==3) {
          A = sqrtm( Sig(0.4,d) )
          B = sqrtm( adj*Sig(0.4,d) )
          X = t( A%*%matrix(rchisq(d*m, df=3), d) )
          Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(mu, d, n) )
          a = bd.test(X, Y, method = "limit")
      }
  }
  # GT
  if (type==5) {
      if (dat==1) {
          A = Sig(0.4,d)
          X = rmvnorm(m, mean = rep(0, d), sigma = A)
          Y = rmvnorm(n, mean = rep(mu, d), sigma = adj*A)
          E = mstree(dist(rbind(X,Y)), 5)
          a = g.tests(E, 1:n, (n+1):(m+n), test.type="g", perm=10000)
      }
      if (dat==2) {
          A = Sig(0.4,d)
          B = Sig(0.4,d)
          X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
          Y = rmvt(n, sigma = adj*A, delta = rep(mu, d), df = 20)
          E = mstree(dist(rbind(X,Y)), 5)
          a = g.tests(E, 1:n, (n+1):(m+n), test.type="g", perm=10000)
      }
      if (dat==3) {
          A = sqrtm( Sig(0.4,d) )
          B = sqrtm( adj*Sig(0.4,d) )
          X = t( A%*%matrix(rchisq(d*m, df=3), d) )
          Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(mu, d, n) )
          E = mstree(dist(rbind(X,Y)), 5)
          a = g.tests(E, 1:n, (n+1):(m+n), test.type="g", perm=10000)
      }
  }
  # CT
  if (type==6) {
      if (dat==1) {
          A = Sig(0.4,d)
          X = rmvnorm(2*m, mean = rep(0, d), sigma = A)
          Y = rmvnorm(2*n, mean = rep(mu, d), sigma = adj*A)
          a = classifier_test(X, Y, split = 0.5)
      }
      if (dat==2) {
          A = Sig(0.4,d)
          B = Sig(0.4,d)
          X = rmvt(2*m, sigma = A, delta = rep(0, d), df = 20)
          Y = rmvt(2*n, sigma = adj*A, delta = rep(mu, d), df = 20)
          a = classifier_test(X, Y, split = 0.5)
      }
      if (dat==3) {
          A = sqrtm( Sig(0.4,d) )
          B = sqrtm( adj*Sig(0.4,d) )
          X = t( A%*%matrix(rchisq(d*2*m, df=3), d) )
          Y = t( B%*%matrix(rchisq(d*2*n, df=3), d) + matrix(mu, d, 2*n) )
          a = classifier_test(X, Y, split = 0.5)
      }
  }
  
  return(a)
}
