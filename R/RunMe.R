install.packages('mvtnorm') 
install.packages('PearsonDS') # for MMD-Pearson
install.packages('ade4')
install.packages('Ball')
install.packages('gTests')
install.packages('Ecume')
install.packages('expm')
install.packages('stats')
install.packages('MASS')
install.packages('abind')

source('kerTests.R')
source('useful_functions.R')
source('simul.R')
source('deep.R')
source('HDMMD.R')
source('AG.R')
source('simul1.R')
source('simul_power.R')
source('simul_power1.R')

###############################################################################
###############################################################################
# Simulations for Section 4 ###################################################
###############################################################################
###############################################################################

#############################################################
############################### Obtain p-value ###############################
#############################################################

# Run: simul(type, m, n, d, dat, a, sigma2, permut)
# type: 1=New tests, 2=MMD Bootstrap, 3=MMD Pearson, 4=BT, 5=GT, 6=CT
# m, n: number of samples
# d: dimension of data
# dat: 1=Gaussian data used in a paper, 2=Multi-t used in a paper, 3=Chi-square used in a paper
# a: mean difference
# sigma2: variance difference
# permut: number of permutation or bootstrap replicates

# In the paper, we use the following parameters:
# Gaussian: m=n=50, d=c(50,100,500,1000), a=c(0.16,0.15,0.1,0.09), sigma2=c(1.11,1.09,1.05,1.04)
# Multi-t: m=n=50, d=c(50,100,500,1000), a=c(0.14,0.15,0.17,0.17), sigma2=c(1.16,1.2,1.2,1.25)
# Chi-square: m=n=50, d=c(50,100,500,1000), a=c(0.29,0.29,0.24,0.24), sigma2=c(1.12,1.11,1.06,1.06)

# Example1: New tests for Gaussian data with changes in both the mean and variance
res = simul(type=1, m=50, n=50, d=500, dat=1, a=0.05, sigma2=1.04, permut=1000)
res$pval$fGPK_appr # p-value of fGPK based on approximation
res$pval$fGPK_perm # p-value of fGPK based on permutation
res$pval$fGPKM_appr # p-value of fGPKM based on approximation
res$pval$fGPKM_perm # p-value of fGPKM based on permutation
res$pval$GPK_perm # p-value of GPK based on permutation

# Example2: MMD-Bootstrap for Gaussian data with changes in both the mean and variance
res = simul(type=2, m=50, n=50, d=500, dat=1, a=0.05, sigma2=1.04, permut=1000)
res # p-value of MMD-Bootstrap

# Example3: MMD-Pearson for Gaussian data with changes in both the mean and variance
res = simul(type=3, m=50, n=50, d=500, dat=1, a=0.05, sigma2=1.04, permut=1000)
res # p-value of MMD-Pearson

# Example4: BT for Gaussian data with changes in both the mean and variance
res = simul(type=4, m=50, n=50, d=500, dat=1, a=0.05, sigma2=1.04, permut=1000)
res$p.value # p-value of BT

# Example5: GT for Gaussian data with changes in both the mean and variance
res = simul(type=5, m=50, n=50, d=500, dat=1, a=0.05, sigma2=1.04, permut=1000)
res$generalized$pval.appr # p-value of GT

# Example6: CT for Gaussian data with changes in both the mean and variance
res = simul(type=6, m=50, n=50, d=500, dat=1, a=0.05, sigma2=1.04, permut=1000)
res$p.value # p-value of CT

#############################################################
############################### Obtain power  #################################
#############################################################

# Run: simul_power(type, m, n, d, dat, a, sigma2, permut, iter, alpha)
# type: 1=New tests, 2=MMD Bootstrap, 3=MMD Pearson, 4=BT, 5=GT, 6=CT
# m, n: number of samples
# d: dimension of data
# dat: 1=Gaussian data used in a paper, 2=Multi-t used in a paper, 3=Chi-square used in a paper
# a: mean difference
# sigma2: variance difference
# permut: number of permutation or bootstrap replicates
# iter: iteration to obtain the empirical power
# alpha: significance level

# In the paper, we use the following parameters:
# Gaussian: m=n=50, d=c(50,100,500,1000), a=c(0.16,0.15,0.1,0.09), sigma2=c(1.11,1.09,1.05,1.04), iter=1000, alpha=0.05
# Multi-t: m=n=50, d=c(50,100,500,1000), a=c(0.14,0.15,0.17,0.17), sigma2=c(1.16,1.2,1.2,1.25),  iter=1000, alpha=0.05
# Chi-square: m=n=50, d=c(50,100,500,1000), a=c(0.29,0.29,0.24,0.24), sigma2=c(1.12,1.11,1.06,1.06),  iter=1000, alpha=0.05

# Example
res = simul_power(type=1, m=50, n=50, d=500, dat=1, a=0.05, sigma2=1.04, permut=1000, iter=1000, alpha=0.05)
res # power of the new tests (%)

###############################################################################
###############################################################################
# Simulations for Section 1.2 and 1.3 #########################################
###############################################################################
###############################################################################

#############################################################
############################### Obtain p-value ###############################
#############################################################

# Run: simul1(type, m, n, d, a, b, permut)
# type: 1=New tests, 2=MMD Bootstrap, 3=AG, 4=ND1, 5=ND2, 6=TM
# m, n: number of samples
# d: dimension of data
# a: mean difference
# b: variance difference
# permut: number of permutation or bootstrap replicates

# In the paper, we use the following parameters:
# Setting1: m=n=50, d=50, a=0.21, b=1
# Setting2: m=n=50, d=50, a=0.21, b=1.04
# Setting3: m=n=50, d=50, a=0, b=1.1

# Example1: New tests for Setting 1
res = simul1(type=1, m=50, n=50, d=50, a=0.21, b=1, permut=1000)
res$pval$fGPK_appr # p-value of fGPK based on approximation
res$pval$fGPK_perm # p-value of fGPK based on permutation
res$pval$fGPKM_appr # p-value of fGPKM based on approximation
res$pval$fGPKM_perm # p-value of fGPKM based on permutation
res$pval$GPK_perm # p-value of GPK based on permutation

# Example2: MMD Bootstrap for Setting 1
res = simul1(type=2, m=50, n=50, d=50, a=0.21, b=1, permut=1000)
res # p-value of MMD Bootstrap

# Example3: AG for Setting 1
res = simul1(type=3, m=50, n=50, d=50, a=0.21, b=1, permut=1000)
res # p-value of AG

# Example4: ND1 for Setting 1
res = simul1(type=4, m=50, n=50, d=50, a=0.21, b=1, permut=1000)
res # p-value of ND1

# Example5: ND2 for Setting 1
res = simul1(type=5, m=50, n=50, d=50, a=0.21, b=1, permut=1000)
res # p-value of ND2

# Example6: TM for Setting 1
res = simul1(type=6, m=50, n=50, d=50, a=0.21, b=1, permut=1000)
res # p-value of TM

#############################################################
############################### Obtain power  #################################
#############################################################

# Run: simul_power1(type, m, n, d, a, b, permut, iter, alpha)
# type: 1=New tests, 2=MMD Bootstrap, 3=AG, 4=ND1, 5=ND2, 6=TM
# m, n: number of samples
# d: dimension of data
# a: mean difference
# b: variance difference
# permut: number of permutation or bootstrap replicates
# iter: iteration to obtain the empirical power
# alpha: significance level

# In the paper, we use the following parameters:
# Setting1: m=n=50, d=50, a=0.21, b=1, iter=1000, alpha=0.05
# Setting2: m=n=50, d=50, a=0.21, b=1.04, iter=1000, alpha=0.05
# Setting3: m=n=50, d=50, a=0, b=1.1, iter=1000, alpha=0.05

# Example
res = simul_power1(type=1, m=50, n=50, d=50, a=0.21, b=1, permut=1000, iter=1000, alpha=0.05)
res # power of the new tests (%)
