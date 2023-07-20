library(mvtnorm)
library(ade4)
library(Ball)
library(gTests)
library(Ecume)
library(expm)

simul = function(type, m, n, d, dat, a, sigma2, permut) {
  # new test
  if (type==1) {
    # Gaussian
    if (dat==1) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A) 
      Y = rmvnorm(n, mean = rep(a, d), sigma = sigma2*A)
      sigma = med_sigma(X,Y)
      res = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
    }
    # t
    if (dat==2) {
      A = Sig(0.4,d)
      B = Sig(0.4,d)
      X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
      Y = rmvt(n, sigma = sigma2*A, delta = rep(a, d), df = 20)
      sigma = med_sigma(X,Y)
      res = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
    }
    # chi
    if (dat==3) {
      A = sqrtm( Sig(0.4,d) )
      B = sqrtm( sigma2*Sig(0.4,d) )
      X = t( A%*%matrix(rchisq(d*m, df=3), d) )
      Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(a, d, n) )
      sigma = med_sigma(X,Y)
      res = kertests(X, Y, sigma, r1=1.2, r2=0.8, perm=permut)
    }
  }
  # MMD Bootstrap
  if (type==2) {
    if (dat==1) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = sigma2*A)
      sigma = med_sigma(X,Y)
      res = mmd(X, Y, sigma, B=permut)
    }
    if (dat==2) {
      A = Sig(0.4,d)
      B = Sig(0.4,d)
      X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
      Y = rmvt(n, sigma = sigma2*A, delta = rep(a, d), df = 20)
      sigma = med_sigma(X,Y)
      res = mmd(X, Y, sigma, B=permut)
    }
    if (dat==3) {
      A = sqrtm( Sig(0.4,d) )
      B = sqrtm( sigma2*Sig(0.4,d) )
      X = t( A%*%matrix(rchisq(d*m, df=3), d) )
      Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(a, d, n) )
      sigma = med_sigma(X,Y)
      res = mmd(X, Y, sigma, B=permut)
    }
  }
  # MMD Pearson
  if (type==3) {
    if (dat==1) {
      A = Sig(0.4,d)
      X = rmvnorm(m, mean = rep(0, d), sigma = A)
      Y = rmvnorm(n, mean = rep(a, d), sigma = sigma2*A)
      sigma = med_sigma(X,Y)
      res = mPearson(X, Y, sigma)
    }
    if (dat==2) {
      A = Sig(0.4,d)
      B = Sig(0.4,d)
      X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
      Y = rmvt(n, sigma = sigma2*A, delta = rep(a, d), df = 20)
      sigma = med_sigma(X,Y)
      res = mPearson(X, Y, sigma)
    }
    if (dat==3) {
      A = sqrtm( Sig(0.4,d) )
      B = sqrtm( sigma2*Sig(0.4,d) )
      X = t( A%*%matrix(rchisq(d*m, df=3), d) )
      Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(a, d, n) )
      sigma = med_sigma(X,Y)
      res = mPearson(X, Y, sigma)
    }
  }
  # BT
  if (type==4) {
      if (dat==1) {
          A = Sig(0.4,d)
          X = rmvnorm(m, mean = rep(0, d), sigma = A)
          Y = rmvnorm(n, mean = rep(a, d), sigma = sigma2*A)
          res = bd.test(X, Y, method = "limit")
      }
      if (dat==2) {
          A = Sig(0.4,d)
          B = Sig(0.4,d)
          X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
          Y = rmvt(n, sigma = sigma2*A, delta = rep(a, d), df = 20)
          res = bd.test(X, Y, method = "limit")
      }
      if (dat==3) {
          A = sqrtm( Sig(0.4,d) )
          B = sqrtm( sigma2*Sig(0.4,d) )
          X = t( A%*%matrix(rchisq(d*m, df=3), d) )
          Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(a, d, n) )
          res = bd.test(X, Y, method = "limit")
      }
  }
  # GT
  if (type==5) {
      if (dat==1) {
          A = Sig(0.4,d)
          X = rmvnorm(m, mean = rep(0, d), sigma = A)
          Y = rmvnorm(n, mean = rep(a, d), sigma = sigma2*A)
          E = mstree(dist(rbind(X,Y)), 5)
          res = g.tests(E, 1:n, (n+1):(m+n), test.type="g")
      }
      if (dat==2) {
          A = Sig(0.4,d)
          B = Sig(0.4,d)
          X = rmvt(m, sigma = A, delta = rep(0, d), df = 20)
          Y = rmvt(n, sigma = sigma2*A, delta = rep(a, d), df = 20)
          E = mstree(dist(rbind(X,Y)), 5)
          res = g.tests(E, 1:n, (n+1):(m+n), test.type="g")
      }
      if (dat==3) {
          A = sqrtm( Sig(0.4,d) )
          B = sqrtm( sigma2*Sig(0.4,d) )
          X = t( A%*%matrix(rchisq(d*m, df=3), d) )
          Y = t( B%*%matrix(rchisq(d*n, df=3), d) + matrix(a, d, n) )
          E = mstree(dist(rbind(X,Y)), 5)
          res = g.tests(E, 1:n, (n+1):(m+n), test.type="g")
      }
  }
  # CT
  if (type==6) {
      if (dat==1) {
          A = Sig(0.4,d)
          X = rmvnorm(2*m, mean = rep(0, d), sigma = A)
          Y = rmvnorm(2*n, mean = rep(a, d), sigma = sigma2*A)
          res = classifier_test(X, Y, split = 0.5)
      }
      if (dat==2) {
          A = Sig(0.4,d)
          B = Sig(0.4,d)
          X = rmvt(2*m, sigma = A, delta = rep(0, d), df = 20)
          Y = rmvt(2*n, sigma = sigma2*A, delta = rep(a, d), df = 20)
          res = classifier_test(X, Y, split = 0.5)
      }
      if (dat==3) {
          A = sqrtm( Sig(0.4,d) )
          B = sqrtm( sigma2*Sig(0.4,d) )
          X = t( A%*%matrix(rchisq(d*2*m, df=3), d) )
          Y = t( B%*%matrix(rchisq(d*2*n, df=3), d) + matrix(a, d, 2*n) )
          res = classifier_test(X, Y, split = 0.5)
      }
  }
  
  return(res)
}
