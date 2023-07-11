# kerTests

This repository provides the _R_ functions for the paper: 

**Song, H.**, Chen, H.    \
  [Generalized kernel two-sample tests.](https://arxiv.org/pdf/2011.06127.pdf) \
  _arXiv_ &ensp; [R package: [kerTests](https://cran.r-project.org/web/packages/kerTests/index.html)]

* The main function, **omni**, provides the omnibus p-value by combining popular normalization strategies, such as none, rarefaction, TSS, CSS, and CLR, based on the Cauchy combination test (Liu et al., 2020).

* The main function utilizes association testing methods: linear regression, ZINQ proposed by Ling et al. (2021), and QRank proposed by Song et al. (2017), which are have been shown to consistently protect type I error and do not depend on the particular choice of normalization.
