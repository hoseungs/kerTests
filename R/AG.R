nnld_MMD2_ts<-function(d2mat,s,permutation,n)#two-sample unbiased est of MMD^2
{
  gkmat<-exp(-d2mat*s)
  
  diag(gkmat)<-0
  N<-nrow(gkmat)
  m<-N-n
  
  v1<-permutation[1:n]
  v2<-permutation[(n+1):N]
  
  MMD2_ts<-sum(gkmat[v1,v1])/(n*(n-1))-2*mean(gkmat[v1,v2])+sum(gkmat[v2,v2])/(m*(m-1));
  return(MMD2_ts)
}

Var_each<-function(gkmat)
{
  n<-dim(gkmat)[1]
  gkmat2<-gkmat^2
  
  diag(gkmat)<-0
  diag(gkmat2)<-0
  
  Comp1<-sum(gkmat2)
  Comp2<-sum((apply(gkmat,MARGIN=2,FUN=sum))^2)-Comp1
  Comp3<-(sum(gkmat))^2-4*Comp2-2*Comp1
  
  Var_each<-Comp1/(n*(n-1))-2*Comp2/(n*(n-1)*(n-2))+Comp3/(n*(n-1)*(n-2)*(n-3))
  
  return(Var_each)
}

Var_ts<-function(d2mat,s)#est of the variance of the test statistic
{
  gkmat<-exp(-d2mat*s)
  return(Var_each(gkmat))
}

mHSIC<-function(gkarray)#unbiased for k=2
{
  # H<-dim(gkarray)[3]
  H<-2
  n<-dim(gkarray)[1]
  for (h in 1:H)
  {
    diag(gkarray[,,h])<-0
  }
  comp1<-sum(apply(gkarray,c(1,2),prod))
  comp2<-sum(apply(apply(gkarray,c(1,3),sum),1,prod))-comp1
  comp3<-prod(apply(gkarray,3,sum))-4*comp2-2*comp1
  
  result<-comp1/(n*(n-1))-2*comp2/(n*(n-1)*(n-2))+comp3/(n*(n-1)*(n-2)*(n-3))
  return(result)
}

mHSIC_alt<-function(gkarray)#slightly biased
{
  H<-dim(gkarray)[3]
  n<-dim(gkarray)[1]
  for (h in 1:H)
  {
    diag(gkarray[,,h])<-0
  }
  comp1<-sum(apply(gkarray,c(1,2),prod))
  comp2<-sum(apply(apply(gkarray,c(1,3),sum),1,prod))
  comp3<-prod(apply(gkarray,3,sum))
  
  result<-comp1/(n*(n-1))-2*comp2/(n*(n-1)^H)+comp3/(n*(n-1))^H
  return(result)
}

nnld_ind<-function(d2array,s,permutation)
{
  gkarray<-exp(-d2array*s)
  
  H<-dim(d2array)[3]
  for (h in 2:H)
  {
    pmt<-permutation[,h-1];
    gkarray[,,h]<-gkarray[pmt,pmt,h]
  }
  
  nnld_ind<-mHSIC(gkarray)
  # nnld_ind<-mHSIC_alt(gkarray)
  
  return(nnld_ind)
}

Var_joint<-function(gkarray)
{
  Var_joint<-prod(apply(gkarray,3,FUN=Var_each));
  return(Var_joint)
}

Var_joint_alt<-function(gkarray)
{
  gkmat<-apply(gkarray,c(1,2),FUN=prod);
  n<-dim(gkmat)[1];
  gkmat2<-gkmat^2;
  
  diag(gkmat)<-0;
  diag(gkmat2)<-0;
  
  Comp1<-sum(gkmat2);
  Comp2<-sum((apply(gkmat,MARGIN=2,FUN=sum))^2)-Comp1;
  Comp3<-(sum(gkmat))^2-4*Comp2-2*Comp1;
  
  Var_joint_alt<-Comp1/(n*(n-1))-2*Comp2/(n*(n-1)*(n-2))+Comp3/(n*(n-1)*(n-2)*(n-3));
  return(Var_joint_alt)
}


