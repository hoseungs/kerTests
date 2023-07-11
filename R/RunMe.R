install.packages('mvtnorm') 
install.packages('PearsonDS') # for MMD-Pearson
install.packages('ade4')
install.packages('Ball')
install.packages('gTests')
install.packages('Ecume')
install.packages('expm')

source('kerTests.R')
source('useful_functions.R')
source('simul.R')



# Run: simul(type, m, n, d, dat, mu, adj, permut)
# type: 1=new tests, 2=MMD Bootstrap, 3=MMD Pearson, 4=BT, 5=GT, 6=CT
# m, n: number of samples
# d: dimension of data
# dat: 1=Gaussian data used in a paper, 2=Log-normal used in a paper, 1=Multi-t used in a paper
# mu: mean difference
# adj: variance difference
# permut: number of permutation or bootstrap replicates



# Example1: new tests for Gaussian data with changes in both the mean and variance
a = simul(type=1, m=100, n=100, d=500, dat=1, mu=0.05, adj=1.04, permut=1000)
a$pval$fGPK_appr # p-value of fGPK based on approximation
a$pval$fGPK_perm # p-value of fGPK based on permutation
a$pval$fGPKM_appr # p-value of fGPKM based on approximation
a$pval$fGPKM_perm # p-value of fGPKM based on permutation
a$pval$GPK_perm # p-value of GPK based on permutation

# Example2: MMD-Bootstrap for Gaussian data with changes in both the mean and variance
a = simul(type=2, m=100, n=100, d=500, dat=1, mu=0.05, adj=1.04, permut=1000)
a # p-value of MMD-Bootstrap


# Example3: MMD-Pearson for Gaussian data with changes in both the mean and variance
a = simul(type=3, m=100, n=100, d=500, dat=1, mu=0.05, adj=1.04, permut=1000)
a # p-value of MMD-Pearson









